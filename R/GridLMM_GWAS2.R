prepMM = function(formula,data,weights = NULL,other_formulas = NULL,
                  relmat = NULL, X = NULL, X_ID = 'ID',proximal_markers = NULL,
                  verbose = TRUE) {
  
  # ensure data has rownames
  if(is.null(rownames(data))) rownames(data) = 1:nrow(data)
  
  # check that data has appropriate columns
  extra_terms = unique(c(X_ID,do.call(c,lapply(other_formulas,function(x) all.vars(x)))))
  if(!all(extra_terms %in% colnames(data))) stop(sprintf("Terms '%s' missing from data",paste(setdiff(extra_terms,colnames(data)),collapse=', ')))
  
  # check that RE's have approporiate levels in data
  RE_levels = list() # a list of levels for each of the random effects
  if(is.null(relmat)) relmat = list()
  for(re in names(relmat)) {
    # check that K is a matrix, then convert to Matrix
    if(is.list(relmat[[re]])){
      if(is.data.frame(relmat[[re]]$K)) relmat[[re]]$K = as.matrix(relmat[[re]]$K)
      if(is.matrix(relmat[[re]]$K)) relmat[[re]]$K = Matrix(relmat[[re]]$K,sparse=T)
      if(is.null(rownames(relmat[[re]]$K))) stop(sprintf('K %s must have rownames',re))
      RE_levels[[re]] = rownames(relmat[[re]]$K)
    } else{
      if(is.data.frame(relmat[[re]])) relmat[[re]] = as.matrix(relmat[[re]])
      if(is.matrix(relmat[[re]])) relmat[[re]] = Matrix(relmat[[re]],sparse=T)
      if(is.null(rownames(relmat[[re]]))) stop(sprintf('K %s must have rownames',re))
      RE_levels[[re]] = rownames(relmat[[re]])
    }
  }
  for(re in names(RE_levels)){
    if(!re %in% colnames(data)) stop(sprintf('Column "%s" required in data',re))
    data[[re]] = as.factor(data[[re]]) # ensure 'data[[re]]' is a factor
    if(!all(data[[re]] %in% RE_levels[[re]])) stop(sprintf('Levels of random effect %s missing.',re))
    data[[re]] = factor(data[[re]],levels = RE_levels[[re]]) # add levels to data[[re]]
  }
  if(!is.null(X_ID)) {
    if(!X_ID %in% colnames(data)) stop(sprintf('X_ID column %s not in data',X_ID))
    data[[X_ID]] = as.factor(data[[X_ID]])
    # if(is.null(rownames(X)) || !all(data[[X_ID]] %in% rownames(X))) stop('X must have rownames and all levels of data[[X_ID]] must be in rownames(X)')
  }
  
  # Use lme4 to evaluate formula in data
  lmod <- lme4::lFormula(formula,data=data,weights=weights,
                         control = lme4::lmerControl(check.nobs.vs.nlev = 'ignore',check.nobs.vs.nRE = 'ignore'))
  
  # Add in other variables from data
  for(term in extra_terms) {
    if(term %in% colnames(lmod$fr) == F) lmod$fr[[term]] = data[match(rownames(lmod$fr),rownames(data)),term]
  }
  
  # compute RE_setup
  RE_terms = lmod$reTrms
  
  if(!is.null(X_ID) && !X_ID %in% names(RE_terms$cnms)) warning(sprintf("No covariance given SNPs specified in error formula. To specify, add a term like (1|%s)",X_ID))
  
  # construct the RE_setup list
  # contains:
  # Z: n x r design matrix
  # K: r x r PSD covariance matrix
  RE_setup = list()
  for(i in 1:length(RE_terms$cnms)){
    term = names(RE_terms$cnms)[i]
    n_factors = length(RE_terms$cnms[[i]])  # number of factors for this grouping factor
    
    # extract combined Z matrix
    combined_Zt = RE_terms$Ztlist[[i]]
    Zs_term = tapply(1:nrow(combined_Zt),gl(n_factors,1,nrow(combined_Zt),labels = RE_terms$cnms[[i]]),function(x) Matrix::t(combined_Zt[x,,drop=FALSE]))
    
    # extract K from relmat. If missing, assign to NULL
    K = NULL
    p = NULL
    p_test = NULL
    if(!is.null(X_ID) && term == X_ID){
      if(term %in% names(relmat)){
        if(verbose) print('using provided RRM.')
        if(is.list(relmat[[term]])){
          if(verbose && !is.null(proximal_markers)) print('Note: X should be centered and scaled as it was for calculating RRM')
          K = relmat[[term]]$K
          if(!is.null(relmat[[term]]$p)){
            p = relmat[[term]]$p
          }
        } else{
          K = relmat[[term]]
        }
      } else{
        if(verbose) print('making RRM matrix')
        p = ncol(X)
        K = tcrossprod(X)/p
        relmat[[term]] = K
      }
      if(!is.null(proximal_markers)) {
        p_test = mean(lengths(proximal_markers))
      }
    } else if(term %in% names(relmat)){
      K = relmat[[term]]
    } 
    
    if(!is.null(K)) {
      if(!all(colnames(Zs_term[[1]]) %in% rownames(K))) stop('rownames of K not lining up with Z')
      K = K[colnames(Zs_term[[1]]),colnames(Zs_term[[1]])]
    } else {
      K = Diagonal(ncol(Zs_term[[1]]),1)
      rownames(K) = colnames(K) = colnames(Zs_term[[1]])
    }
    
    
    # make an entry in RE_setup for each random effect
    for(j in 1:n_factors){
      # name of variance component
      name = term
      if(n_factors > 1) name = paste(name,RE_terms$cnms[[i]][[j]],sep='.')
      while(name %in% names(RE_setup)) name = paste0(name,'.1') # hack for when same RE used multiple times
      
      # Z matrix
      Z = as(Zs_term[[j]],'dgCMatrix')
      
      
      RE_setup[[name]] = list(
        term = term,
        Z = Z,
        K = K,
        p = p,
        p_test = p_test
      )
    }
    
    # add names to RE_setup if needed
    n_RE = length(RE_setup)
    for(i in 1:n_RE){
      if(is.null(names(RE_setup)[i]) || names(RE_setup)[i] == ''){
        names(RE_setup)[i] = paste0('RE.',i)
      }
    }
    
  }
  
  return(list(lmod = lmod, RE_setup = RE_setup))
}

setup_Grid2 = function(RE_names,h2_step,h2_start = NULL){
  
  # -------- Form matrix of h2s to test ---------- #
  n_RE = length(RE_names)
  
  if(length(h2_step) < n_RE){
    if(length(h2_step) != 1) stop('Must provide either 1 h2_step parameter, or 1 for each random effect')
    h2_step = rep(h2_step,n_RE)
  }
  if(is.null(names(h2_step))) {
    names(h2_step) = RE_names
  }
  
  if(length(h2_start) != length(RE_names)) stop("Wrong length of h2_start")
  if(is.null(h2_start)) h2_start = rep(0,length(RE_names))
  if(is.null(names(h2_start))) names(h2_start) = RE_names
  
  h2s_matrix = expand.grid(lapply(RE_names,function(re) h2_start[[re]] + seq(-1,2,by = h2_step[[re]])))
  colnames(h2s_matrix) = RE_names
  h2s_matrix = h2s_matrix[rowSums(h2s_matrix<0) == 0,,drop=FALSE]
  h2s_matrix = t(h2s_matrix[rowSums(h2s_matrix) < 1,,drop=FALSE])
  colnames(h2s_matrix) = NULL
  
  return(h2s_matrix)
}

GridLMM_GWAS2 = function(formula,test_formula,reduced_formula,data,weights = NULL,
                         X,X_ID = 'ID',
                         centerX = FALSE,scaleX = FALSE,fillNAX = FALSE,X_map = NULL, relmat = NULL,
                         h2_step = 0.01, h2_start = NULL, h2_start_tolerance = 0.001,max_steps = 100, method = c('REML','ML','BF'), algorithm = c('Fast','Full'),
                         inv_prior_X = NULL,target_prob = 0.99,
                         downdate_REs = NULL, proximal_markers = NULL, proximal_Xs = NULL,
                         V_setup = NULL, save_V_folder = NULL, 
                         diagonalize=T,
                         mc.cores = my_detectCores(),clusterType = 'mclapply',
                         verbose=T) {
  
  # Evaluate argument options
  # -------- Evaluate argument options ---------- #
  method = match.arg(method)
  algorithm = match.arg(algorithm)
  
  # -------- Re-scale SNPs ---------- #
  X = scale_SNPs(X,centerX,scaleX,fillNAX)
  
  # -------- prep Mixed Models ---------- #
  MM = prepMM(formula,data,weights,other_formulas = list(test_formula,reduced_formula),
              relmat,X,X_ID,proximal_markers,verbose)
  lmod = MM$lmod
  RE_setup = MM$RE_setup
  
  Y = matrix(lmod$fr[,1])
  colnames(Y) = 'y'
  X_cov = lmod$X
  data = lmod$fr

  # -------- Prepare tests ---------- #
  # -------- Test and reduced models ---------- #
  # test model should is applied to each SNP. The intercept is the main effect, and any additional terms are multiplied by the SNP
  mm_test = model.matrix(test_formula,droplevels(data))
  mm_reduced = model.matrix(reduced_formula,droplevels(data))
  
    
  # -------- prep Markers and proximal_Xs ---------- #
  # check that proximal_Xs is appropriate
  downdate_REs = NULL
  if(!is.null(proximal_markers)) {
    if(length(RE_setup) > 1) {
      # downdate_REs is a vector of length == length(RE_setup) with 1s identifying REs that will be downdated
      downdate_REs = sapply(names(RE_setup),function(x) strsplit(x,'.',fixed=TRUE)[[1]]) == X_ID
      downdate_REs = as.numeric(downdate_REs)
      
      if(is.null(proximal_Xs)) {
        # use X_list_full
        if(ncol(mm_test) != length(downdate_REs)) stop("Can't tell what markers to downdate from each RE")
        proximal_Xs = downdate_REs
        proximal_Xs[proximal_Xs == 1] = 1:ncol(mm_test)
      } else if(is.list(proximal_Xs)) {
        if(length(proximal_Xs) != sum(downdate_REs) && length(proximal_Xs) != length(downdate_REs)) stop("Wrong length of proximal_Xs")
        if(!is.null(names(proximal_Xs))) {
          if(!all(names(proximal_Xs) %in% names(RE_setup)[downdate_REs == 1])) {
            stop(sprintf("proximal_Xs must be named: %s",paste(names(RE_setup),collapse=',')))
          }
          proximal_Xs = lapply(names(RE_setup),function(x) proximal_Xs[[x]])
        }
      }
    } else{
      if(is.null(proximal_Xs)) {
        proximal_Xs = 1
      } else {
        if(length(proximal_Xs) > 1) stop("Wrong length of proximal_Xs")
      }
    }
  }
  
  # model.matrix used to expand X to size of data based on X_ID
  Z_X = model.matrix(formula(sprintf("~0+%s",X_ID)),droplevels(data))
  colnames(Z_X) = sub(X_ID,'',colnames(Z_X),fixed = T)
  if(is.null(rownames(X))) stop(sprintf('X must have rownames that correspond with column %s in data',X_ID))
  stopifnot(all(colnames(Z_X) %in% rownames(X)))
  
  X = Z_X %*% X[colnames(Z_X),]
  if(is.list(proximal_Xs)) {
    for(i in 1:length(proximal_Xs)) {
      if(is.null(rownames(proximal_Xs[[i]]))) stop(sprintf('proximal_Xs must have rownames that correspond with column %s in data',X_ID))
      proximal_Xs[[i]] = Z_X %*% proximal_Xs[[i]][colnames(Z_X),]
    }
  }
  
  # -------- prep V_setup ---------- #
  if(is.null(V_setup)) {
    V_setup = make_V_setup(RE_setup,weights,diagonalize,svd_K = TRUE,drop0_tol = 1e-10,save_V_folder,verbose)
  } 
  
  # -------- setup prior for X ------ #
  if(method == 'BF') {
    if(is.null(inv_prior_X)) {
      warning('Calculating Bayes Factors with impropper prior on X might be unreliable')
      inv_prior_X = rep(0,ncol(X_cov) + ncol(mm_test))
    } else if(length(inv_prior_X) == 1) {
      inv_prior_X = c(rep(0,ncol(X_cov)),rep(inv_prior_X,ncol(mm_test)))
    } else if(length(inv_prior_X) == ncol(mm_test)) {
      inv_prior_X = c(rep(0,ncol(X_cov)),inv_prior_X)
    } else if(length(inv_prior_X) != (ncol(X_cov)+ncol(mm_test))) {
      stop('Wrong length of inv_prior_X')
    }
  } else{
    inv_prior_X = rep(0,ncol(X_cov) + ncol(mm_test))
  }
  
  # -------- get starting value by EMMAX ---------- #
  # evaluate likelihoods on a grid, going until all LL's stop increasing
  if(!is.null(Y) && is.null(h2_start)) {# && method != 'BF') {
    if(verbose) print('Estimating h2_start via null model')
    if(ncol(Y) > 1) stop('EMMAX_start only implemented for a single response')
    if(method == 'BF') {
      null_Bayes = GridLMM_posterior(formula,data,V_setup = V_setup,h2_divisions=3,mc.cores = mc.cores,verbose = verbose)
      h2_start = matrix(colSums(null_Bayes$h2s_results$posterior*null_Bayes$h2s_results[,1:length(V_setup$ZKZts),drop=FALSE]),nr=1)
    } else{
      null_ML = GridLMM_ML(formula,data,V_setup = V_setup,tolerance = h2_start_tolerance,mc.cores = mc.cores,verbose=verbose)
      ML = REML = FALSE
      if(method == 'ML') ML = TRUE
      if(method == 'REML' || method == 'BF') REML = TRUE
      h2_start = get_current_h2s(null_ML$results,names(V_setup$ZKZts),ML = ML,REML=REML)
    }
    if(verbose) {
      print('GridLMM_posterior h2s:')
      print(h2_start)
    }
  }
  
  if(algorithm == 'Full') {
    h2_start = t(setup_Grid2(names(RE_setup),h2_step,h2_start))
  }
  
  
  setup = list(
    Y = Y,
    X_cov = X_cov,
    X = X,
    h2_start = h2_start,
    V_setup = V_setup,
    mm_test = mm_test,
    mm_reduced = mm_reduced,
    proximal_markers = proximal_markers,
    proximal_Xs = proximal_Xs
  )
  
  
  results = run_GWAS_fast2(Y,X_cov,X,mm_test,mm_reduced,inv_prior_X,V_setup,h2_start,h2_step,target_prob,proximal_markers, proximal_Xs, method,mc.cores,verbose)
  
  return(list(
    results = results,
    setup = setup
    ))
}


run_GWAS_fast2 = function(Y,X_cov,X,mm_test,mm_reduced = NULL,inv_prior_X = NULL,V_setup,h2_start=NULL,h2_step,target_prob = 0.99,proximal_markers=NULL, proximal_Xs=NULL, 
               method = c('REML','ML','BF'), mc.cores=my_detectCores(),verbose = FALSE)
{
  
  # -------- Evaluate argument options ---------- #
  method = match.arg(method)
  
  # -------- setup prior for X ------ #
  if(method == 'BF') {
    if(is.null(inv_prior_X)) {
      warning('Calculating Bayes Factors with impropper prior on X might be unreliable')
      inv_prior_X = rep(0,ncol(X_cov) + ncol(mm_test))
    } else if(length(inv_prior_X) == 1) {
      inv_prior_X = c(rep(0,ncol(X_cov)),rep(inv_prior_X,ncol(mm_test)))
    } else if(length(inv_prior_X) == ncol(mm_test)) {
      inv_prior_X = c(rep(0,ncol(X_cov)),inv_prior_X)
    } else if(length(inv_prior_X) != (ncol(X_cov)+ncol(mm_test))) {
      stop('Wrong length of inv_prior_X')
    }
  } else{
    inv_prior_X = rep(0,ncol(X_cov) + ncol(mm_test))
  }
  
  X_list_full = lapply(1:ncol(mm_test),function(j) mm_test[,j] * X)
  if(ncol(mm_reduced) > 0 && method %in% c('ML','BF')) {
    X_list_reduced = lapply(1:ncol(mm_reduced),function(j) mm_reduced[,j] * X)
  } else{
    X_list_reduced = NULL
  }
  
  # -------- rotate data ---------- #
  # speeds to everything if a single random effect

  Qt = V_setup$Qt
  if(!is.null(Qt)){
    Y <- as.matrix(Qt %*% Y)
    X_cov <- as.matrix(Qt %*% X_cov)
    test_names = colnames(X_list_full[[1]])
    X_list_full <- premultiply_list_of_matrices(Qt, X_list_full)
    colnames(X_list_full[[1]]) = test_names
    if(length(X_list_reduced)>0) {
      reduced_names = colnames(X_list_reduced[[1]])
      X_list_reduced <- premultiply_list_of_matrices(Qt, X_list_reduced)
      colnames(X_list_reduced[[1]]) = reduced_names
    }
    if(is.list(proximal_Xs)) {
      proximal_Xs <- premultiply_list_of_matrices(Qt,proximal_Xs)
    }
  }
  
  # -------- prep_h2s ---------- #
  RE_names = names(V_setup$ZKZts)
  if(is.matrix(h2_start)){
    if(!ncol(h2_start) == length(RE_names)) stop('wrong length of h2_start provided')
    if(is.null(colnames(h2_start))) colnames(h2_start) = RE_names
  } else if(!length(h2_start) == length(RE_names)) {
    stop('wrong length of h2_start provided')
    if(is.null(names(h2_start))) names(h2_start) = RE_names
  }
  
  # -------- run GWAS_fast ---------- #
  # evaluate likelihoods on a grid, going until all LL's stop increasing
  
  results = fit_GWAS_fast2(Y,X_cov,X_list_full,h2_start,h2_step,V_setup,inv_prior_X,target_prob,proximal_markers,proximal_Xs,method,verbose,mc.cores)
  
  # Do tests
  if(method == 'REML') {
    F_hat_cols = colnames(results)[grep('F.',colnames(results),fixed=T)]
    p_value_REML = do.call(cbind,lapply(F_hat_cols,function(col) p_value = pf(results[,col],1,nrow(Y) - results$Df_X - ncol(X_cov),lower.tail=F)))
    if(length(F_hat_cols) > 1) {
      colnames(p_value_REML) = sprintf('p_value_REML.%d',1:length(F_hat_cols))
    } else{
      colnames(p_value_REML) = 'p_value_REML'
    }
    results = data.frame(results,p_value_REML)
  } else {
    if(is.null(proximal_Xs) || !is.list(proximal_Xs)) {
      if(!is.list(proximal_Xs)) {
        proximal_Xs_new = list()
        proximal_Xs_new[proximal_Xs] = X_list_full
        proximal_Xs = proximal_Xs_new
        rm(X_list_full)
        gc()
      }
    }
    if(method == 'ML') {
      results_reduced = fit_GWAS_fast2(Y,X_cov,X_list_reduced,h2_start,h2_step,V_setup,inv_prior_X,target_prob,proximal_markers,proximal_Xs,method,verbose)
      for(trait in unique(results$Trait)){
        index_full = results$Trait == trait
        index_reduced = results_reduced$Trait == trait
        results$ML_Reduced_logLik[index_full] = results_reduced$ML_logLik[index_reduced]
        results$Reduced_Df_X[index_full] = results_reduced$Df_X[index_reduced]
        results$Reduced_X_id[index_full] = results_reduced$X_id[index_reduced]
        
        results$p_value_ML = with(results,pchisq(2*(ML_logLik - ML_Reduced_logLik),Df_X - Reduced_Df_X,lower.tail = F))
      }
    }
    if(method == 'BF'){
      results_reduced = fit_GWAS_fast2(Y,X_cov,X_list_reduced,h2_start,h2_step,V_setup,inv_prior_X,target_prob,proximal_markers,proximal_Xs,method,verbose)
      for(trait in unique(results$Trait)){
        index_full = results$Trait == trait
        index_reduced = results_reduced$Trait == trait
        results$log_posterior_factor_reduced = results_reduced$log_posterior_factor[index_reduced]
        results$BF[index_full] = results$log_posterior_factor - results$log_posterior_factor_reduced # Note: need to account for different V_beta
        results$Reduced_Df_X[index_full] = results_reduced$Df_X[index_reduced]
        if(!is.null(inv_prior_X)){
          if(length(inv_prior_X) == 1) {
            results$BF[index_full] = results$BF[index_full] - log(inv_prior_X)/2
          } else{
            results$BF[index_full] = results$BF[index_full] - sum(log(tail(inv_prior_X,n = results$Df_X[1] - results$Reduced_Df_X[1])))/2
          }
        }
      }
    }
  }
  results$ML_index = c()
  results$REML_index = c()
  
  return(results)
}


fit_GWAS_fast2 = function(Y,X_cov,X_list,h2_start,h2_step,V_setup,inv_prior_X = NULL,target_prob = 0.99,
                          proximal_markers = NULL,proximal_Xs = NULL,
                          method = 'REML',verbose=F,mc.cores = 1) {
  # calculates LL's for each test over a grid of h2s 
  # starts at h2_start, and then moves out in equidistant circles. 
  # Only tests where the LL increases in one circle are re-tested in the next.
  # continues until the grid is fully covered, or no tests increase in LL
  Y = as.matrix(Y)
  storage.mode(Y) = 'double'
  if(is.null(colnames(Y))) colnames(Y) = 1:ncol(Y)
  
  if(method == 'ML'){
    ML = TRUE
    REML = FALSE
    BF = FALSE
  }
  if(method == 'REML'){
    ML = FALSE
    REML = TRUE
    BF = FALSE
  }
  if(method == 'BF'){
    ML = FALSE
    REML = FALSE
    BF = TRUE
  }
  downdate_Xs = NULL
  
  null_model = FALSE
  if(is.null(X_list) || length(X_list) == 0) {
    if(is.null(proximal_markers)) null_model = TRUE
  } else if(is.null(colnames(X_list[[1]])) & ncol(X_list[[1]])>0) colnames(X_list[[1]]) = 1:ncol(X_list[[1]])
  
  if(is.null(inv_prior_X)){
    inv_prior_X = c(rep(0,ncol(X_cov)),rep(0,length(X_list)))
  }
  
  RE_names = colnames(h2_start)
  if(is.null(RE_names)) RE_names = names(h2_start)
  if(null_model) {
    active = TRUE
    h2s_to_test = as.matrix(h2_start)
    tested_h2s = c()
    results = c()
    
    n_steps = 1
    while(active && nrow(h2s_to_test) > 0) {   
      if(verbose) print(sprintf('step: %d, num_h2s: %d, num active: %d',n_steps,nrow(h2s_to_test), sum(active)))
      inner_cores = 1
      outer_cores = mc.cores
      registerDoParallel(outer_cores)
      results_list = foreach(h2s = iter(h2s_to_test,by='row')) %dopar% {
        h2s = h2s[1,]
        chol_V_setup = make_chol_V_setup(V_setup,h2s)
        chol_V = chol_V_setup$chol_V
        V_log_det <- chol_V_setup$V_log_det
        calc_LL_parallel2(Y,X_cov,X_list,h2s,chol_V,V_log_det,inv_prior_X[1:ncol(X_cov)],
                          NULL,V_setup$n_SNPs_downdated_RRM,REML,BF,inner_cores,c())
      }
      tested_h2s = rbind(tested_h2s,h2s_to_test)
      if(length(results) == 0) {
        results = compile_results(results_list)
      } else{
        results = compile_results(c(list(results),results_list))
      }
      
      if(n_steps > 1) {
        active = FALSE
        if(ML && !is.null(results$ML_index) && results$ML_index > 1) active = TRUE
        if(REML && !is.null(results$REML_index) && results$REML_index > 1) active = TRUE
        if(BF && results$BF_new_posterior_mass > 1-target_prob) active = TRUE
      }
      if(active) {
        results$n_steps = n_steps
      } 
      
      if(!BF) {
        current_h2s = get_current_h2s(results,RE_names,ML,REML)
      } else{
        current_h2s = tested_h2s
      }
      h2s_to_test = get_h2s_to_test(current_h2s,tested_h2s,h2_step,ML,REML)
      if(ncol(h2s_to_test) == 0) break
      n_steps = n_steps + 1   
    }
    
    results$Df_X = rep(0,ncol(Y))
    
  } else{
    
    h2s_to_test = as.matrix(h2_start)
    tested_h2s = c()
    
    if(is.null(proximal_markers)) {
      active_tests = rep(TRUE,ncol(X_list[[1]]))
    } else{
      active_tests = rep(TRUE,length(proximal_markers))
    }
    n_steps = 1
    results = c()
    
    while(sum(active_tests) > 0 && nrow(h2s_to_test) > 0) {
      if(verbose) print(sprintf('step: %d, num_h2s: %d, num active: %d',n_steps,nrow(h2s_to_test), sum(active_tests)))
      active_X_list = which(active_tests)
      
      if(!is.null(proximal_markers)) {
        if(is.list(proximal_Xs)){
          downdate_Xs = build_downdate_Xs(1:length(proximal_Xs),proximal_Xs,proximal_markers[active_X_list])
        } else{
          downdate_Xs = build_downdate_Xs(proximal_Xs,X_list,proximal_markers[active_X_list])
        }
      }
      
      # test each h2 in h2s_to_test
      inner_cores = 1
      outer_cores = max(1,min(mc.cores,nrow(h2s_to_test)))
      if(outer_cores == 1) inner_cores = mc.cores
      registerDoParallel(outer_cores)
      results_list = foreach(h2s = iter(h2s_to_test,by='row')) %dopar% {
        h2s = h2s[1,]
        chol_V_setup = make_chol_V_setup(V_setup,h2s)
        chol_V = chol_V_setup$chol_V
        V_log_det <- chol_V_setup$V_log_det
        calc_LL_parallel2(Y,X_cov,X_list,h2s,chol_V,V_log_det,inv_prior_X,
                          downdate_Xs,V_setup$n_SNPs_downdated_RRM,REML,BF,inner_cores,active_X_list)
      }
      tested_h2s = rbind(tested_h2s,h2s_to_test)
      
      # for each test, identify the best h2s value and associated statistics
      if(length(results) == 0) {
        results = compile_results(results_list)
      } else{
        results = compile_results(c(list(results),results_list))
      }
      
      # determine which tests haven't maxed out yet
      if(n_steps > 1) {
        active_tests[] = FALSE
        if(ML && !is.null(results$ML_index)) active_tests[!is.na(results$ML_index) & results$ML_index > 1] = TRUE
        if(REML && !is.null(results$REML_index)) active_tests[!is.na(results$REML_index) & results$REML_index > 1] = TRUE
        if(BF) active_tests[results$BF_new_posterior_mass > 1-target_prob] = TRUE
      }
      results$n_steps[active_tests] = n_steps
      
      if(!BF) {
        current_h2s = get_current_h2s(results,RE_names,ML,REML)
      } else{
        current_h2s = tested_h2s
      }
      h2s_to_test = get_h2s_to_test(current_h2s,tested_h2s,h2_step,ML,REML)
      if(nrow(h2s_to_test) == 0) break
      n_steps = n_steps + 1
      
    }
    results$Df_X = rep(length(X_list),ncol(Y))
    if(length(X_list) == 0) results$Df_X = 0
    
  }
  return(results)
}




calc_LL2 = function(Y,X_cov,X_list,h2s,chol_Vi,V_log_det,inv_prior_X,
                    downdate_Xs = NULL,n_SNPs_downdated_RRM = NULL,REML = TRUE, BF = TRUE,active_X_list = NULL){  
  if(is.null(active_X_list)) {
    if(!is.null(downdate_Xs)) {
      active_X_list = 1:length(downdate_Xs)
    } else if(length(X_list) == 0) {
      active_X_list = integer()
    } else {
      active_X_list = 1:ncol(X_list[[1]])
    }
  }
  m = ncol(Y)
  n = nrow(Y)
  if(is.null(downdate_Xs)) {
    SSs <- GridLMM_SS_matrix(Y,chol_Vi,X_cov,X_list,active_X_list,inv_prior_X);
  } else{
    downdate_weights = h2s/unlist(n_SNPs_downdated_RRM)
    if(length(downdate_weights) != length(downdate_Xs[[1]]$downdate_Xi)) stop("Wrong length of downdate weights")
    chol_Vi = as.matrix(chol_Vi)
    SSs <- GridLMM_SS_downdate_matrix(Y,chol_Vi,X_cov,X_list,active_X_list,downdate_Xs,downdate_weights,inv_prior_X);
  }
  log_LL = get_LL(SSs,X_cov,X_list,active_X_list,n,m,ML=TRUE,REML = REML,BF=BF)
  
  # collect results in a table
  # beta_hats and F_hats are lists of matrices, with coefficients as columns and traits as rows. 
  # We want to organize into a single matrix with coefficients as columns and each trait sequentially, with all tests for each trait consecutive
  
  # p = max(1,max(ncol(X_list[[1]]),length(downdate_weights_i)))*m
  # ML_h2 = as.matrix(t(h2s))[rep(1,p),,drop=FALSE]
  n_tests = nrow(log_LL$ML)
  if(is.null(n_tests)) n_tests = 1
  p = n_tests*m
  ML_h2 = as.matrix(t(h2s))[rep(1,p),,drop=FALSE]
  colnames(ML_h2) = paste(colnames(ML_h2),'ML',sep='.')
  
  results_i = data.frame(Trait = rep(colnames(Y),each = n_tests),
                         X_ID = NA,
                         ML_logLik = c(log_LL$ML),
                         ML_h2,
                         stringsAsFactors = F)
  if(length(active_X_list) == n_tests) {
    if(length(X_list) > 0) {
      results_i$X_ID = colnames(X_list[[1]])[active_X_list]
    } else {
      results_i$X_ID = names(downdate_Xs)[active_X_list]
    }
  } else{
    results_i$X_ID = 'NULL'
  }
  
  b_cov = ncol(X_cov)
  b_x = length(X_list)
  if(length(X_list) == 0) b_x = 0
  b = b_cov + b_x
  # if(ncol(X_list[[1]]) == n_tests) {
  # beta_hats = do.call(rbind,lapply(1:m,function(i) t(log_LL$beta_hat[(i-1)*b + b_cov + 1:b_x,,drop=FALSE])))
  beta_hats = do.call(rbind,lapply(1:m,function(i) t(log_LL$beta_hat[(i-1)*b + 1:b,,drop=FALSE])))
  colnames(beta_hats) = paste('beta',1:ncol(beta_hats),sep='.')
  results_i = data.frame(results_i,beta_hats)
  # }
  
  if(REML) {
    REML_h2 = as.matrix(t(h2s))[rep(1,p),,drop=FALSE]
    colnames(REML_h2) = paste(colnames(REML_h2),'REML',sep='.')
    results_i = data.frame(results_i,REML_logLik = c(log_LL$REML),REML_h2)
    if(b_x > 0) {
      F_hats = do.call(rbind,lapply(1:m,function(i) t(log_LL$F_hat[(i-1)*b + b_cov + 1:b_x,,drop=FALSE])))
      colnames(F_hats) = paste('F',1:ncol(F_hats),sep='.')
      results_i = data.frame(results_i,F_hats)
    }
  }
  if(BF){
    results_i$log_posterior_factor = c(t(log_LL$log_posterior_factor))
    results_i$RSSs = c(t(log_LL$RSSs))
  }
  return(results_i)
}

calc_LL_parallel2 = function(Y,X_cov,X_list,h2s,chol_V,V_log_det,inv_prior_X,
                             downdate_Xs = NULL,n_SNPs_downdated_RRM = NULL,REML = TRUE,BF = FALSE,
                             mc.cores = 1,active_X_list = NULL) {
  
  if(mc.cores == 1 || length(X_list) == 0) {
    return(calc_LL2(Y,X_cov,X_list,h2s,chol_V,V_log_det,inv_prior_X,downdate_Xs,n_SNPs_downdated_RRM,REML,BF,active_X_list))
  }
  # make X_list_active, divide into chunks
  X_list_active = NULL
  total_tests = 0
  if(!is.null(X_list)) {
    X_list_active = lapply(X_list,function(x) x[,active_X_list,drop=FALSE])
    total_tests = ncol(X_list_active[[1]])
    if(!is.null(downdate_Xs) && (total_tests != length(downdate_Xs))) stop("Wrong length of downdate_Xs")
  } else{
    total_tests = length(downdate_Xs)
  }
  # X_list_names = colnames(X_list[[1]])
  registerDoParallel(mc.cores)
  chunkSize = total_tests/mc.cores
  chunks = 1:total_tests
  chunks = split(chunks, ceiling(seq_along(chunks)/chunkSize))
  X_list_sets = lapply(chunks,function(i) lapply(X_list_active,function(Xi) Xi[,i,drop=FALSE]))
  if(!is.null(downdate_Xs)) {
    downdate_Xs_sets = lapply(chunks,function(i) downdate_Xs[i])
    results = foreach(X_list_i = iter(X_list_sets),downdate_Xs_i = iter(downdate_Xs_sets),.combine = 'rbind') %dopar% {
      calc_LL2(Y,X_cov,X_list_i,h2s,chol_V,V_log_det,inv_prior_X,downdate_Xs_i,n_SNPs_downdated_RRM,REML,BF)
    }
  } else{
    results = foreach(X_list_i = iter(X_list_sets),.combine = 'rbind') %dopar% {
      calc_LL2(Y,X_cov,X_list_i,h2s,chol_V,V_log_det,inv_prior_X,downdate_Xs,n_SNPs_downdated_RRM,REML,BF)
    }
  }
  return(results)
}


