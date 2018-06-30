

make_V_list = function(RE_setup,
                       weights = NULL,
                       h2_divisions,   
                       h2_EMMAX = NULL,
                       save_V_list = NULL, # character vector giving folder name to save V_list
                       diagonalize = TRUE,
                       svd_K = TRUE, # should the svd of one of the random effect covariances be calculated?
                       drop0_tol = 1e-10, # tolerance for setting a value to zero in a sparse matrix
                       mc.cores = parallel::detectCores(),
                       clusterType = 'mclapply',
                       verbose = T
){
  # either returns a folder containing 1) each chol_V_inv file (containing a chol_V_inv, h2_index, and V_log_det) and a chol_V_inv_setup file, 2) the h2s_matrix, and 3) the Qt matrix
  # or, a list containing these three items in memory
  
  if(clusterType == 'SNOW') {
    cl = makeCluster(mc.cores)
  } else{
    cl = mc.cores
  }
  
  # ------------------------------------ #
  # ---- Check RE_setup ---------------- #
  # ------------------------------------ #
  
  # ensure both Z and K are provided for each random effect
  for(re in seq_along(RE_setup)){
    RE_setup[[re]] = within(RE_setup[[re]],{
      if(!'Z' %in% ls()){
        Z = diag(1,nrow(K))
      }
      if(!'K' %in% ls()){
        K = diag(1,ncol(Z))
      }
    })
  }
  
  
  # ------------------------------------ #
  # ---- Check weights ----------------- #
  # ------------------------------------ #
  
  n = nrow(RE_setup[[1]]$Z)
  if(is.null(weights)) weights = rep(1,n)
  if(length(weights) != n) stop('weights must have same length as data')
  
  # ------------------------------------ #
  # --------------- h2s ---------------- #
  # ------------------------------------ #
  
  # table of possible h2s for each random effect
  #   These are percentages of the total residual variance accounted for by each random effect
  #   Each column is a set of percentages, the sum of which must be less than 1 (so that Ve is > 0)
  # can specify different levels of granularity for each random effect
  
  RE_names = names(RE_setup)
  n_RE = length(RE_names)
  
  if(length(h2_divisions) < n_RE){
    if(length(h2_divisions) != 1) stop('Must provide either 1 h2_divisions parameter, or 1 for each random effect')
    h2_divisions = rep(h2_divisions,n_RE)
  }
  if(is.null(names(h2_divisions))) {
    names(h2_divisions) = RE_names
  }
  h2s_matrix = expand.grid(lapply(RE_names,function(re) seq(0,1,length = h2_divisions[[re]]+1)))
  colnames(h2s_matrix) = RE_names
  h2s_matrix = t(h2s_matrix[rowSums(h2s_matrix) < 1,,drop=FALSE])
  if(!is.null(h2_EMMAX)) {
    if(length(h2_EMMAX) != nrow(h2s_matrix)) stop('wrong length h2_EMMAX vector provided')
    if(!is.null(names(h2_EMMAX))){
      if(!all(rownames(h2s_matrix) %in% names(h2_EMMAX))) stop('missing h2s in h2_EMMAX')
      h2s_matrix = cbind(h2s_matrix,h2_EMMAX[rownames(h2s_matrix)])
    } else{
      h2s_matrix = cbind(h2s_matrix,h2_EMMAX)
    }
  }
  colnames(h2s_matrix) = NULL
  
  # ------------------------------------ #
  # ----Precalculate Sigmas ------------ #
  # ------------------------------------ #
  
  # do calculations in several chunks
  group_size = 2*mc.cores
  n_groups = ceiling(ncol(h2s_matrix)/group_size)
  col_groups = tapply(1:ncol(h2s_matrix),gl(n_groups,group_size,ncol(h2s_matrix)),function(x) x)
  
  if(verbose) {
    print(sprintf("Pre-calculating random effect inverse matrices for %d sets of random effect weights", ncol(h2s_matrix)))
    pb = txtProgressBar(min=0,max = length(col_groups),style=3)
  }
  
  if(n_RE == 1 && diagonalize == TRUE && all(weights == 1)){
    # only useful if n_RE == 1 and there are no weights
    Z = RE_setup[[1]]$Z
    if(svd_K == TRUE && nrow(RE_setup[[1]]$K) < nrow(Z)) {
      # a faster way of taking the SVD of ZKZt, particularly if ncol(Z) < nrow(Z). Probably no benefit if ncol(K) > nrow(Z)
      svd_K1 = svd(RE_setup[[1]]$K)
      qr_ZU = qr(Z %*% svd_K1$u)
      R_ZU = drop0(qr.R(qr_ZU,complete=F),tol=drop0_tol)
      Q_ZU = drop0(qr.Q(qr_ZU,complete=T),tol=drop0_tol)
      RKRt = R_ZU %*% diag(svd_K1$d) %*% t(R_ZU)
      svd_RKRt = svd(RKRt)
      RKRt_U = svd_RKRt$u
      if(ncol(Q_ZU) > ncol(RKRt_U)) RKRt_U = bdiag(RKRt_U,diag(1,ncol(Q_ZU)-ncol(RKRt_U)))
      Qt = t(Q_ZU %*% RKRt_U)
    } else{
      ZKZt = Z %*% RE_setup[[1]]$K %*% t(Z)
      result = svd(ZKZt)
      Qt = t(result$u)
    }
    Qt = as(drop0(as(Qt,'dgCMatrix'),tol = drop0_tol),'dgCMatrix')
    QtZ_matrices = lapply(RE_setup,function(re) Qt %*% re$Z)
  } else{
    Qt = NULL
    QtZ_matrices = lapply(RE_setup,function(re) re$Z)
  }
  QtZ = do.call(cbind,QtZ_matrices[RE_names])
  QtZ = as(QtZ,'dgCMatrix')
  
  ZKZts = list()
  for(i in 1:n_RE){
    ZKZts[[i]] = as(forceSymmetric(drop0(QtZ_matrices[[i]] %*% RE_setup[[i]]$K %*% t(QtZ_matrices[[i]]),tol = drop0_tol)),'dgCMatrix')
  }
  
  chol_V_list = c()
  
  if(is.character(save_V_list)){
    try(dir.create(save_V_list))
    system(sprintf('rm %s/*',save_V_list),ignore.stderr = T)
  }
  
  # in preparation for down-dating, multiply each ZKZt by p/(p-pj), where p is the number of SNPs that went into the RRM, and pj is the (mean) number of SNPs per down-date operation
  n_SNPs_downdated_RRM = lapply(RE_setup,function(x) {
    if(!is.null(x$p) && !is.null(x$p_test)) return(x$p - x$p_test)
    return(0)
  })
  downdate_ratios = c()
  for(i in 1:n_RE){
    if(is.null(RE_setup[[i]]$p)) {
      downdate_ratios[i] = 1
    } else if(is.null(RE_setup[[i]]$p_test)){
      downdate_ratios[i] = RE_setup[[i]]$p/(RE_setup[[i]]$p - 1) # assume p_test == 1 if not provided
    } else{
      downdate_ratios[i] = RE_setup[[i]]$p/(RE_setup[[i]]$p - RE_setup[[i]]$p_test) 
    }
  }
  
  V_list_setup = list(
    h2s_matrix = h2s_matrix,
    Qt = Qt,
    n_SNPs_downdated_RRM = n_SNPs_downdated_RRM
  )
  
  registerDoParallel(cl)
  for(group in col_groups){
    chol_V_list = c(chol_V_list,foreach::foreach(h2_index=group) %dopar% {
      h2s = h2s_matrix[,h2_index]
      V = diag((1-sum(h2s))/weights)  # include weights here. This should be the only place they are needed.
      for(i in 1:n_RE) V = V + h2s[i] * ZKZts[[i]] * downdate_ratios[i]
      chol_V = chol(V)
      if(length(chol_V@i) > prod(dim(chol_V))/4) chol_V = as.matrix(chol_V) # no reason to store as sparse matrix if not sparse
      V_log_det = 2*sum(log(diag(chol_V)))
      return(list(h2_index = h2_index,chol_V = chol_V, V_log_det = V_log_det))
    })
    if(is.character(save_V_list)){
      for(i in 1:length(chol_V_list)){
        chol_V = chol_V_list[[i]]
        chol_V_file = sprintf('%s/chol_V_%05d.rds',save_V_list,group[i])
        saveRDS(chol_V,file = chol_V_file)
        V_list_setup$chol_V_list = c(V_list_setup$chol_V_list,chol_V_file)
      }
      chol_V_list = c()
    }
    if(verbose) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
  }
  if(verbose) close(pb)
  
  if(is.character(save_V_list)){
    setup_file = sprintf('%s/V_list_setup.rds',save_V_list)
    saveRDS(V_list_setup,file = setup_file)
    V_list_setup = list(
      setup_file = setup_file,
      chol_V_files = V_list_setup$chol_V_list
    )
  } else{
    V_list_setup$chol_V_list = chol_V_list
  }
  return(V_list_setup)
}

GridLMM_test <- function(
  Y, # (n x m) data vector
  X_cov = NULL, # (n x p) matrix of covariates constant for each model
  X_list_full, # list of model matrices to use as the full model
  X_list_reduced = NULL, # list of model matrices to use as a reduced model. The covariates should be removed from both X_list_full and X_list_reduced. Can be a single matrix if the reduced model is the same for all tests
  V_list_setup, # either a list returned by \link{make_V_list}, or a file name pointing to the \code{V_list_setup} object
  inv_prior_X = NULL,
  downdate_Xs = NULL, # list of matrices for each random effect and each test to downdate covariance matrices for each test. Optional. Expected if V_listchol_V_list calculated with down-date weights
  REML = TRUE,
  BF = TRUE,
  mc.cores = parallel::detectCores(),
  clusterType = 'mclapply'
) {
  Y = as.matrix(Y)
  X_cov = as.matrix(X_cov)
  n = nrow(Y)
  m = ncol(Y)
  if(is.null(X_cov)) X_cov = matrix(1,n,1)
  
  if(is.null(X_list_reduced)){
    X_list_reduced = list(matrix(0,nrow=n,ncol=0))
  }
  
  # check that number of columns of X_list_full tests is constant
  if(length(X_list_full) > 1 && !var(sapply(X_list_full,ncol)) == 0) stop('all tests in X_list_full must have same number of columns')
  
  # check that X_list_reduced is compatible with X_list_full (either length==0 or each element has same number of columns as X_list_full[[1]])
  if(ncol(X_list_reduced[[1]]) > 0){
    if(any(sapply(X_list_reduced,ncol) != ncol(X_list_full[[1]]))) stop('number of tests in X_list_full is different than in X_list_reduced')
  }
  
  # check that downdate_Xs is compatible with X_list_full
  if(!is.null(downdate_Xs)){
    if(any(sapply(downdate_Xs,function(Xs) length(Xs)) != ncol(X_list_full[[1]]))){
      stop('number of downdate_Xs must equal number of tests (ie length(X_list_full)')
    }
  }
  
  test_names = colnames(X_list_full[[1]])
  if(!is.null(downdate_Xs)) test_names = names(downdate_Xs[[1]])
  
  if(is.character(V_list_setup)) {
    V_list_setup = readRDS(V_list_setup)
  } else if('setup_file' %in% names(V_list_setup)) {
    V_list_setup = readRDS(V_list_setup$setup)
  }
  Qt = V_list_setup$Qt
  
  # do this here because premultiply_list_of_matrices does a deep copy of downdate_Xs  - I don't want it done twice!
  if(!is.null(Qt)){
    Y <- as.matrix(Qt %*% Y)
    X_cov <- as.matrix(Qt %*% X_cov)
    X_list_full <- premultiply_list_of_matrices(Qt, X_list_full)
    if(ncol(X_list_full[[1]]) > 0) colnames(X_list_full[[1]]) = test_names
    X_list_reduced <- premultiply_list_of_matrices(Qt, X_list_reduced)
    if(ncol(X_list_reduced[[1]]) > 0) colnames(X_list_reduced[[1]]) = test_names
    if(!is.null(downdate_Xs)) {
      downdate_Xs <- lapply(downdate_Xs,function(Xi) premultiply_list_of_matrices(Qt, Xi))
      names(downdate_Xs[[1]]) = test_names
    }
  }
  
  logLL_full <- GridLMM_logLL(
    Y = Y,
    X_cov = X_cov,
    X_list = X_list_full,
    V_list_setup = V_list_setup,
    inv_prior_X = inv_prior_X,
    REML = REML,
    BF = BF,
    downdate_Xs,
    mc.cores = mc.cores,
    clusterType = clusterType
  )
  
  logLL_reduced <- GridLMM_logLL(
    Y = Y,
    X_cov = X_cov,
    X_list = X_list_reduced,
    V_list_setup = V_list_setup,
    inv_prior_X = inv_prior_X,
    REML = FALSE,
    BF = BF,
    downdate_Xs,
    mc.cores = mc.cores,
    clusterType = clusterType
  )
  
  results = logLL_full
  for(trait in unique(results$Trait)){
    index_full = results$Trait == trait
    index_reduced = logLL_reduced$Trait == trait
    results$ML_Reduced_logLik[index_full] = logLL_reduced$ML_logLik[index_reduced]
    results$Reduced_Df_X[index_full] = logLL_reduced$Df_X[index_reduced]
    results$Reduced_X_ID[index_full] = logLL_reduced$X_ID[index_reduced]
    results$Reduced_ML_index = logLL_reduced$ML_index[index_reduced]
    
    results$p_value_ML = with(results,pchisq(2*(ML_logLik - ML_Reduced_logLik),Df_X - Reduced_Df_X,lower.tail = F))
  }
  if(REML) {
    F_hat_cols = colnames(results)[grep('F.',colnames(results),fixed=T)]
    p_value_REML = do.call(cbind,lapply(F_hat_cols,function(col) p_value = pf(results[,col],1,nrow(Y) - results$Df_X - ncol(X_cov),lower.tail=F)))
    if(length(F_hat_cols) > 1) {
      colnames(p_value_REML) = sprintf('p_value_REML.%d',1:length(F_hat_cols))
    } else{
      colnames(p_value_REML) = 'p_value_REML'
    }
    results = data.frame(results,p_value_REML)
  }
  if(BF){
    for(trait in unique(results$Trait)){
      index_full = results$Trait == trait
      index_reduced = logLL_reduced$Trait == trait
      results$log_posterior_factor_reduced = logLL_reduced$log_posterior_factor[index_reduced]
      results$BF[index_full] = results$log_posterior_factor - results$log_posterior_factor_reduced # Note: need to account for different V_beta
      if(!is.null(inv_prior_X)){
        if(length(inv_prior_X) == 1) {
          results$BF[index_full] = results$BF[index_full] - log(inv_prior_X)/2
        } else{
          results$BF[index_full] = results$BF[index_full] - sum(log(tail(inv_prior_X,n = results$Df_X[index_full] - results$Reduced_Df_X[index_full])))/2
        }
      }
    }
  }
  return(results)
}

GridLMM_logLL <- function(
  Y, # (n x m) data vector
  X_cov, # (n x p) matrix of covariates constant for each model
  X_list, # list of design matrices. Will calculate the logLL for each, maximized over \sigma^2
  V_list_setup, # either a list returned by \link{make_V_list}, or a file name pointing to the \code{V_list_setup} object
  inv_prior_X = NULL, # prior on the precision of each predictor
  REML = TRUE, # should REML likelihood and Wald F be calculated? Not needed for LRT
  BF = TRUE,
  downdate_Xs = NULL, #list of matrices for each random effect and each test to downdate covariance matrices for each test. Optional. Expected if V_list calculated with down-date weights
  mc.cores = parallel::detectCores(),
  clusterType
) {
  if(clusterType == 'SNOW') {
    cl = makeCluster(mc.cores)
  } else{
    cl = mc.cores
  }
  
  Y = as.matrix(Y)
  storage.mode(Y) = 'double'
  m = ncol(Y)
  if(is.null(colnames(Y))) colnames(Y) = 1:m
  if(is.null(colnames(X_list[[1]])) & ncol(X_list[[1]])>0) colnames(X_list[[1]]) = 1:ncol(X_list[[1]])
  test_names = colnames(X_list[[1]])
  if(!is.null(downdate_Xs)) test_names = names(downdate_Xs[[1]])
  
  h2s_matrix             = V_list_setup$h2s_matrix
  n_SNPs_downdated_RRM   = V_list_setup$n_SNPs_downdated_RRM # used to calcualte downdate_weights based on downdate_Xs
  chol_V_list            = V_list_setup$chol_V_list
  
  if(is.null(inv_prior_X) || ncol(X_list[[1]]) == 0){
    inv_prior_X = c(rep(0,ncol(X_cov)),rep(0,length(X_list)))
  } else if(length(inv_prior_X) == 1){
    inv_prior_X = c(rep(0,ncol(X_cov)),rep(inv_prior_X,length(X_list)))
  } else if(length(inv_prior_X) == length(X_list)){
    inv_prior_X = c(rep(0,ncol(X_cov)),inv_prior_X)
  } else {
    stop('wrong length of inv_prior_X')
  }
  
  if(!is.null(downdate_Xs)){
    if(length(downdate_Xs) != nrow(h2s_matrix)) stop('need downdate matrices for each random effect')
  }
  
  registerDoParallel(cl)
  results_list <- foreach::foreach(chol_V_setup = chol_V_list) %dopar% {
    if(is.character(chol_V_setup)) {
      chol_V_setup <- readRDS(chol_V_setup)
    }
    h2s <- h2s_matrix[,chol_V_setup$h2_index]
    chol_V = chol_V_setup$chol_V
    V_log_det <- chol_V_setup$V_log_det
    calc_LL(Y,X_cov,X_list,h2s,chol_V,V_log_det,inv_prior_X,downdate_Xs,n_SNPs_downdated_RRM,REML)
  }
  if(inherits(cl,"cluster")) stopCluster(cl)
  
  results = compile_results(results_list)
  results$Df_X = length(X_list)
  if(ncol(X_list[[1]]) == 0) {
    results$Df_X = 0
  }
  
  return(results)
}  
