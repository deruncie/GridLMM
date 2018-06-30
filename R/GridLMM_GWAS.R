GridLMM_GWAS = function(error_model,test_model,reduced_model,data,Y = NULL, weights = NULL,
                       X,X_ID = 'ID',centerX = TRUE,scaleX = FALSE,fillNAX = FALSE,X_map = NULL, relmat = NULL,inv_prior_X = NULL,
                       h2_divisions,EMMAX_start = TRUE,h2_EMMAX = NULL, REML = TRUE,ML = TRUE,BF = FALSE,proximal_matrix = NULL, 
                       RE_setup = NULL, V_list = NULL, save_V_list = NULL,downdate_Xs = NULL,
                       diagonalize=T,svd_K = T,drop0_tol = 1e-10,mc.cores = my_detectCores(),clusterType = 'mclapply',chunkSize = 10000,verbose=T) {
  
  # -------- check terms in models ---------- #
  terms = c(all.vars(error_model), all.vars(test_model),all.vars(reduced_model))
  if(!all(terms %in% colnames(data))) {
    missing_terms = terms[!terms %in% colnames(data)]
    stop(sprintf('terms %s missing from data',paste(missing_terms,sep=', ')))
  }
  if(!is.null(X_ID)) data[[X_ID]] = factor(data[[X_ID]])
  
  # -------- Response ---------- #
  if(is.null(Y)){
    if(length(error_model) == 3){
      # if(length(model) == 2) then there is no response
      Y = as.matrix(data[,all.vars(error_model[[2]])])
    }
  }
  if(nrow(Y) != nrow(data)) stop('nrow(Y) != nrow(data)')
  # remove NAs in observations. Any observation with any NA in any trait will be dropped
  NA_obs = rowSums(is.na(Y))
  Y = Y[!NA_obs,,drop=FALSE]
  data = data[!NA_obs,,drop=FALSE]
  n = nrow(data)
  
  # -------- Constant Fixed effects ---------- #
  X_cov = model.matrix(nobars(error_model),data)
  linear_combos = caret::findLinearCombos(X_cov)
  if(!is.null(linear_combos$remove)) {
    cat(sprintf('dropping column(s) %s to make covariates full rank\n',paste(linear_combos$remove,sep=',')))
    X_cov = X_cov[,-linear_combos$remove]
  }
  if(any(is.na(X_cov))) stop('Missing values in covariates')
  
  # -------- Re-scale SNPs ---------- #
  X = scale_SNPs(X,centerX,scaleX,fillNAX)
  if(is.null(colnames(X))) colnames(X) = 1:ncol(X)
  
  
  # -------- Random effects ---------- #
  if(is.null(RE_setup)) {
    RE_setup = make_RE_setup(model = error_model,data,relmat = relmat,X = X,X_ID = X_ID,proximal_matrix = proximal_matrix,verbose=verbose)
  }
  
  if(!is.null(Y) && EMMAX_start && !BF) {
    # if reporting BF's, don't use EMMAX_start
    if(!is.null(h2_EMMAX)) stop('h2_EMMAX should be NULL if EMMAX_start == TRUE')
    if(ncol(Y) > 1) stop('EMMAX_start only implemented for a single response')
    # if(!all(names(RE_terms$cnms) %in% names(RE_setup))) stop('lme4qtl failed. Please run estimate h2_EMMAX directly')
    relmats = lapply(names(RE_setup),function(term) {
      K = RE_setup[[term]]$K
      diag(K) = diag(K) + 1e-10
      K
    })  # NOTE: adding small element to diagonal to "fix" eigenvalues < 0
    names(relmats) = names(RE_setup)
    emmax <- lme4qtl::relmatLmer(error_model,data=data,relmat = relmats,weights = weights)
    vars = as.data.frame(lme4::VarCorr(emmax))$vcov
    h2_EMMAX = vars/sum(vars)
    h2_EMMAX = h2_EMMAX[-length(h2_EMMAX)]
    names(h2_EMMAX) = NULL
    if(verbose) print(sprintf('EMMAX h2s: %s',paste(signif(h2_EMMAX,2),collapse=', ')))
  }
  if(BF) h2_EMMAX = NULL
  
  
  
  # -------- make downdate_Xs ---------- #
  if(is.null(downdate_Xs) && !is.null(proximal_matrix)) {
    re_X_ID = names(RE_setup)[which(sapply(RE_setup,function(x) !is.null(x$p)))]
    downdate_X_matrices = lapply(re_X_ID,function(x) as.matrix(RE_setup[[x]]$Z %*% X[colnames(RE_setup[[x]]$Z),]))
    names(downdate_X_matrices) = re_X_ID
    downdate_Xs = lapply(re_X_ID,function(x) lapply(1:ncol(X),function(j){
      SNPs = which(proximal_matrix[j,] != 0)
      downdate_X_matrices[[x]][,SNPs,drop=FALSE]
    }))
    names(downdate_Xs) = re_X_ID
    names(downdate_Xs[[1]]) = colnames(X)
    rm(downdate_X_matrices)
  } 
  
  # -------- make V_list ---------- #
  if(is.null(V_list)) {
    V_list_setup = make_V_list(RE_setup,
                               weights,
                               h2_divisions,  
                               h2_EMMAX,
                               save_V_list,
                               diagonalize,
                               svd_K,
                               drop0_tol,
                               mc.cores,
                               verbose)
  }
  
  
  # -------- Prepare tests ---------- #
  # model.matrix used to expand X to size of data based on X_ID
  Z_X = model.matrix(formula(sprintf("~0+%s",X_ID)),droplevels(data))
  colnames(Z_X) = sub(X_ID,'',colnames(Z_X),fixed = T)
  if(is.null(rownames(X))) stop(sprintf('X must have rownames that correspond with column %s in data',X_ID))
  stopifnot(all(colnames(Z_X) %in% rownames(X)))
  
  # -------- Test and reduced models ---------- #
  # test model should is applied to each SNP. The intercept is the main effect, and any additional terms are multiplied by the SNP
  mm_test = model.matrix(test_model,droplevels(data))
  mm_reduced = model.matrix(reduced_model,droplevels(data))
  
  if(!is.null(Y)) {
    n_tests = ncol(X)
    n_chunks = ceiling(n_tests/chunkSize)
    if(verbose) {
      print("Running GWAS")
      pb = txtProgressBar(min=0,max = n_chunks,style=3)
    }
    results = foreach(Xi = ichunk(iter(X,by='column',drop=FALSE,mode='numeric'),chunkSize = chunkSize),.combine='rbind') %do% {
      Xi = do.call(cbind,Xi)
      ZXi = Z_X %*% Xi[colnames(Z_X),]
      X_list_full_i = lapply(1:ncol(mm_test),function(j) mm_test[,j] * ZXi)
      colnames(X_list_full_i[[1]]) = colnames(Xi)
      if(ncol(mm_reduced) > 0) {
        X_list_reduced_i = lapply(1:ncol(mm_reduced),function(j) mm_reduced[,j] * ZXi)
        colnames(X_list_reduced_i[[1]]) = colnames(Xi)
      } else{
        X_list_reduced_i = NULL
      }
      results_i = GridLMM_test(Y,X_cov,X_list_full_i,X_list_reduced_i,V_list_setup,inv_prior_X,downdate_Xs,REML=REML,BF=BF,mc.cores,clusterType)
      if(verbose) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
      results_i
    }
    if(verbose) close(pb)
  }
  
  setup = list(
    RE_setup = RE_setup,
    proximal_matrix = proximal_matrix,
    X_cov = X_cov,
    Z_X = Z_X,
    mm_test = mm_test,
    mm_reduced = mm_reduced,
    V_list_setup = V_list_setup,
    downdate_Xs = downdate_Xs
  )
  if(is.null(Y)) return(setup)
  
  # add map info
  if(!is.null(X_map)){
    index = match(results$X_ID,X_map$snp)
    if(!is.null(index) && length(index) == nrow(results)) {
      results = data.frame(X_map[index,],results)
    }
  }
  
  return(list(
    results = results,
    setup = setup
  ))
  
}