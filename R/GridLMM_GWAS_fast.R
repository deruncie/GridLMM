GridLMM_GWAS_fast = function(error_model,test_model,reduced_model,data,Y = NULL, weights = NULL,
                            X,X_ID = 'ID',centerX = FALSE,scaleX = FALSE,fillNAX = FALSE,X_map = NULL, relmat = NULL,
                            h2_start = NULL, h2_step = 0.01, h2_start_tolerance = 0.001, max_steps = 100, method = 'ML',
                            inv_prior_X = NULL,target_prob = 0.99,
                            proximal_matrix = NULL, downdate_Xs = NULL,
                            V_setup = NULL, save_V_folder = NULL, 
                            diagonalize=T,svd_K = T,drop0_tol = 1e-10,mc.cores = my_detectCores(),clusterType = 'mclapply',chunkSize = 10000,verbose=T) {
  
  # -------- check terms in models ---------- #
  terms = c(all.vars(error_model), all.vars(test_model),all.vars(reduced_model))
  if(!all(terms %in% colnames(data))) {
    missing_terms = terms[!terms %in% colnames(data)]
    stop(sprintf('terms %s missing from data',paste(missing_terms,sep=', ')))
  }
  if(!is.null(X_ID)) data[[X_ID]] = factor(data[[X_ID]])
  
  # -------- Response ---------- #
  n = nrow(data)
  if(is.null(Y)){
    if(length(error_model) == 3){
      # if(length(model) == 2) then there is no response
      Y = as.matrix(data[,all.vars(error_model[[2]])])
    }
  }
  
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
  
  
  if(is.null(V_setup)) {
    # -------- Random effects ---------- #
    RE_setup = make_RE_setup(model = error_model,data,relmat = relmat,X = X,X_ID = X_ID,proximal_matrix = proximal_matrix,verbose=verbose)
    V_setup = make_V_setup(RE_setup,weights,diagonalize,svd_K,drop0_tol,save_V_folder,verbose)
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
  
  # -------- get starting value by EMMAX ---------- #
  # evaluate likelihoods on a grid, going until all LL's stop increasing
  if(!is.null(Y) && is.null(h2_start)) {# && method != 'BF') {
    if(verbose) print('Estimating h2_start via null model')
    if(ncol(Y) > 1) stop('EMMAX_start only implemented for a single response')
    if(method == 'BF') {
      null_Bayes = GxElmm(error_model,data,V_setup = V_setup,h2_divisions=3,mc.cores = mc.cores)
      h2_start = matrix(colSums(null_Bayes$h2s_results$posterior*null_Bayes$h2s_results[,1:length(V_setup$ZKZts),drop=FALSE]),nr=1)
    } else{
      null_ML = GxElmm_ML(error_model,data,V_setup = V_setup,tolerance = h2_start_tolerance,mc.cores = mc.cores)
      ML = REML = FALSE
      if(method == 'ML') ML = TRUE
      if(method == 'REML' || method == 'BF') REML = TRUE
      h2_start = get_current_h2s(null_ML$results,names(V_setup$ZKZts),ML = ML,REML=REML)
    }
    if(verbose) {
      print('GxElmm h2s:')
      print(h2_start)
    }
  }
  
  # -------- make downdate_Xs ---------- #
  if(is.null(downdate_Xs) && !is.null(proximal_matrix)) {
    RE_setup = V_setup$RE_setup
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
  
  setup = list(
    Y = Y,
    X_cov = X_cov,
    h2_start = h2_start,
    V_setup = V_setup,
    Z_X = Z_X,
    mm_test = mm_test,
    mm_reduced = mm_reduced,
    downdate_Xs = downdate_Xs
  )
  if(is.null(Y)) return(setup)
  
  results = run_GWAS_fast(X,setup,h2_step,max_steps,
                          inv_prior_X,target_prob,
                          proximal_matrix, downdate_Xs,
                          method, mc.cores,verbose)
  
  
  
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

run_GWAS_fast = function(X,setup,h2_step = 0.01,max_steps = 100,
                         inv_prior_X = NULL,target_prob = 0.99,
                         proximal_matrix = NULL, downdate_Xs = NULL,
                         method = 'ML',mc.cores = 1,verbose = TRUE) {
  
  Y <- setup$Y
  X_cov <- setup$X_cov
  h2_start = setup$h2_start
  V_setup = setup$V_setup
  RE_setup = V_setup$RE_setup
  Z_X = setup$Z_X
  mm_test = setup$mm_test
  mm_reduced = setup$mm_reduced
  
  if(is.null(colnames(X))) colnames(X) = 1:ncol(X)
  
  # -------- Prepare tests ---------- #
  ZX = Z_X %*% X[colnames(Z_X),]
  rm(X)  # cleanup to recover memory
  X_list_full = lapply(1:ncol(mm_test),function(j) mm_test[,j] * ZX)
  if(ncol(mm_reduced) > 0 && method %in% c('ML','BF')) {
    X_list_reduced = lapply(1:ncol(mm_reduced),function(j) mm_reduced[,j] * ZX)
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
    if(!is.null(downdate_Xs)) downdate_Xs <- lapply(downdate_Xs,function(Xi) premultiply_list_of_matrices(Qt, Xi))
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
  
  results = fit_GWAS_fast(Y,X_cov,X_list_full,h2_start,h2_step,V_setup,inv_prior_X,target_prob,downdate_Xs,method,verbose,mc.cores)
  
  # Do tests
  if(method == 'ML') {
    results_reduced = fit_GWAS_fast(Y,X_cov,X_list_reduced,h2_start,h2_step,V_setup,inv_prior_X,target_prob,downdate_Xs,method,verbose)
    for(trait in unique(results$Trait)){
      index_full = results$Trait == trait
      index_reduced = results_reduced$Trait == trait
      results$ML_Reduced_logLik[index_full] = results_reduced$ML_logLik[index_reduced]
      results$Reduced_Df_X[index_full] = results_reduced$Df_X[index_reduced]
      results$Reduced_X_id[index_full] = results_reduced$X_id[index_reduced]
      
      results$p_value_ML = with(results,pchisq(2*(ML_logLik - ML_Reduced_logLik),Df_X - Reduced_Df_X,lower.tail = F))
    }
  }
  if(method == 'REML') {
    F_hat_cols = colnames(results)[grep('F.',colnames(results),fixed=T)]
    p_value_REML = do.call(cbind,lapply(F_hat_cols,function(col) p_value = pf(results[,col],1,nrow(Y) - results$Df_X - ncol(X_cov),lower.tail=F)))
    if(length(F_hat_cols) > 1) {
      colnames(p_value_REML) = sprintf('p_value_REML.%d',1:length(F_hat_cols))
    } else{
      colnames(p_value_REML) = 'p_value_REML'
    }
    results = data.frame(results,p_value_REML)
  }
  if(method == 'BF'){
    results_reduced = fit_GWAS_fast(Y,X_cov,X_list_reduced,h2_start,h2_step,V_setup,inv_prior_X,target_prob,downdate_Xs,method,verbose)
    for(trait in unique(results$Trait)){
      index_full = results$Trait == trait
      index_reduced = results_reduced$Trait == trait
      results$log_posterior_factor_reduced = results_reduced$log_posterior_factor[index_reduced]
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
  results$ML_index = c()
  results$REML_index = c()
  
  return(results)
}


fit_GWAS_fast = function(Y,X_cov,X_list,h2_start,h2_step,V_setup,inv_prior_X = NULL,target_prob = 0.99,downdate_Xs=NULL,method = 'ML',verbose=F,mc.cores = 1) {
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
  
  null_model = FALSE
  if(is.null(X_list) || length(X_list) == 0) {
    X_list = list(matrix(0,nrow = nrow(Y),ncol = 0))
    if(is.null(downdate_Xs)) null_model = TRUE
  }
  if(is.null(colnames(X_list[[1]])) & ncol(X_list[[1]])>0) colnames(X_list[[1]]) = 1:ncol(X_list[[1]])
  
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
      print(sprintf('step: %d, num_h2s: %d, num active: %d',n_steps,nrow(h2s_to_test), sum(active)))
      inner_cores = 1
      outer_cores = mc.cores
      registerDoParallel(outer_cores)
      results_list = foreach(h2s = iter(h2s_to_test,by='row')) %dopar% {
        h2s = h2s[1,]
        chol_V_setup = make_chol_V_setup(V_setup,h2s)
        chol_V = chol_V_setup$chol_V
        V_log_det <- chol_V_setup$V_log_det
        calc_LL_parallel(Y,X_cov,X_list,h2s,chol_V,V_log_det,inv_prior_X,downdate_Xs,V_setup$n_SNPs_downdated_RRM,REML,BF,inner_cores)
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
    
    if(is.null(downdate_Xs)) {
      active_tests = rep(TRUE,ncol(X_list[[1]]))
    } else{
      active_tests = rep(TRUE,length(downdate_Xs[[1]]))
    }
    n_steps = 1
    results = c()
    
    while(sum(active_tests) > 0 && nrow(h2s_to_test) > 0) {
      print(sprintf('step: %d, num_h2s: %d, num active: %d',n_steps,nrow(h2s_to_test), sum(active_tests)))
      
      active_X_list = which(active_tests)
      
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
        calc_LL_parallel(Y,X_cov,X_list,h2s,chol_V,V_log_det,inv_prior_X,downdate_Xs,V_setup$n_SNPs_downdated_RRM,REML,BF,inner_cores,active_X_list)
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
    if(ncol(X_list[[1]]) == 0) results$Df_X = 0
    
  }
  return(results)
}

