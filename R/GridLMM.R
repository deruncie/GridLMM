GridLMM_posterior = function(model,data,Y = NULL, X = NULL, weights = NULL,relmat = NULL,  
                  h2_divisions = 10,
                  h2_prior = function(h2s,n) 1/n, a = 0, b = 0, inv_prior_X = 0,
                  target_prob = 0.99, # want this much probability to be distributed over at least thresh_nonzero grid squares 
                  thresh_nonzero = 10, thresh_nonzero_marginal = 0,
                  V_setup = NULL, save_V_folder = NULL, # character vector giving folder name to save V_list
                  diagonalize=T,svd_K = T,drop0_tol = 1e-10,mc.cores = my_detectCores(),verbose=T) {
  
  # -------- check terms in models ---------- #
  terms = c(all.vars(model))
  if(!all(terms %in% colnames(data))) {
    missing_terms = terms[!terms %in% colnames(data)]
    stop(sprintf('terms %s missing from data',paste(missing_terms,sep=', ')))
  }
  
  # -------- Response ---------- #
  n = nrow(data)
  if(is.null(Y)){
    if(length(model) == 3){
      # if(length(model) == 2) then there is no response
      Y = as.matrix(data[,all.vars(model[[2]])])
      storage.mode(Y) = 'numeric'
    }
  }
  
  # -------- Constant Fixed effects ---------- #
  X_cov = model.matrix(nobars(model),data)
  if(!is.null(X)) X_cov = cbind(X_cov,X)
  linear_combos = caret::findLinearCombos(X_cov)
  if(!is.null(linear_combos$remove)) {
    cat(sprintf('dropping column(s) %s to make covariates full rank\n',paste(linear_combos$remove,sep=',')))
    X_cov = X_cov[,-linear_combos$remove]
  }
  if(any(is.na(X_cov))) stop('Missing values in covariates')
  p = ncol(X_cov)
  if(is.null(inv_prior_X)) {
    inv_prior_X = rep(0,p)
  } else if(length(inv_prior_X) == 1) {
    inv_prior_X = rep(inv_prior_X,p)
  } else if(length(inv_prior_X) != p) {
    stop('wrong length of inv_prior_X')
  }
  
  # -------- Random effects ---------- #
  RE_setup = make_RE_setup(model = model,data = data,relmat = relmat)
  
  # -------- prep V ---------- #
  if(is.null(V_setup))  V_setup = make_V_setup(RE_setup,weights,diagonalize,svd_K,drop0_tol,save_V_folder,verbose)
  
  # -------- rotate data ---------- #
  # speeds to everything if a single random effect
  Qt = V_setup$Qt
  if(!is.null(Qt)){
    Y <- as.matrix(Qt %*% Y)
    X_cov <- as.matrix(Qt %*% X_cov)
  }
  
  # -------- prep_h2s ---------- #
  RE_names = names(RE_setup)
  n_RE = length(RE_names)
  
  step_size = 1/h2_divisions
  h2s_to_test = make_h2s_matrix(RE_names,seq(0,1,by=step_size))
  
  # while the number of grid points accounting for target_prob % of the posterior is less than thresh_nonzero, keep dividing the h2s_matrix intervals smaller and smaller around the top regions
  h2s_solutions = c()
  tested_h2s = c()
  num_top_grid = 0
  # recover()
  # h2s_to_test = make_h2s_matrix(RE_names,seq(0,1,by=0.01))
  while(nrow(h2s_to_test)>0) {
    if(verbose) print(sprintf('step_size: %s, non_zero: %d, num_to_test: %d', step_size,num_top_grid, nrow(h2s_to_test)))
    
    # -------- calculate posterior of h2s ----------- #
    registerDoParallel(mc.cores)
    new_h2s_solutions = foreach(h2s = iter(h2s_to_test,by = 'row')) %do% {
      chol_V_setup = make_chol_V_setup(V_setup,unlist(h2s))
      if(inherits(chol_V_setup$chol_V,'Matrix')){ 
        chol_Vinv = t(solve(chol_V_setup$chol_V))
      } else{
        chol_Vinv = t(backsolve(chol_V_setup$chol_V,diag(1,ncol(chol_V_setup$chol_V))))
      }
      V_log_det <- chol_V_setup$V_log_det
      
      SSs <- GridLMM_SS_dense_c(Y,as.matrix(chol_Vinv),X_cov, list(matrix(0,n,0)),integer(),inv_prior_X,V_log_det)
      
      mu_star = SSs$beta_hats
      V_star_inv = matrix(0,length(mu_star),length(mu_star))
      V_star_inv[lower.tri(V_star_inv,diag=T)] = SSs$V_star_L
      V_star_inv = tcrossprod(V_star_inv)
      a_star = a + n/2
      b_star = b + 1/2*SSs$RSSs
      V_star_log_det = -SSs$V_star_inv_log_det
      
      # solve integral over s2, beta
      # constant 1
      # log_c1 = a*log(b) - ((n/2 + p/2)*log(2*pi) + V_log_det/2 + sum(log(inv_prior_X))/2 + lgamma(a))
      # 
      #   # constant 2
      # log_c2 = a_star * log(b_star) - (p/2*log(2*pi) + V_star_log_det/2 + lgamma(a_star))
      # 
      # log_LL = log_c1 - log_c2
      
      # alternate version only with terms that remain in posterior
      log_LL = -V_log_det/2 - a_star * log(b_star) +  V_star_log_det/2
      
      # calculate numerator of posterior
      log_post_num = log_LL # + log(h2_prior(h2s,nrow(h2s_matrix)))  # will add this later
      
      # ML = calc_ML(SSs,n)
      # SSs$log_det_X = log_det_of_XtX(X_cov,list(matrix(0,n,0)))
      # REML = calc_REML(SSs,n,1,1)
      
      return(list(log_post_num = log_post_num,
                  h2s = h2s,
                  mu_star = mu_star,
                  V_star_inv = V_star_inv,
                  a_star = a_star,
                  b_star = b_star,
                  V_log_det = V_log_det,
                  V_star_log_det = V_star_log_det
                  # ,ML = ML,
                  # REML = REML
      ))
    }
    
    h2s_solutions = c(new_h2s_solutions,h2s_solutions)
    
    log_post_num = sapply(h2s_solutions,function(x) x$log_post_num + log(h2_prior(x$h2s,1/step_size)))
    log_post_num[is.na(log_post_num)] = -Inf
    
    max_log_post_num = max(log_post_num)
    norm_factor = max_log_post_num + log(sum(exp(log_post_num - max_log_post_num)))
    h2s_results = data.frame(rbind(h2s_to_test,tested_h2s),posterior = exp(log_post_num - norm_factor))
    tested_h2s = rbind(h2s_to_test,tested_h2s)
    rownames(tested_h2s) = NULL
    
    
    # sort by posterior
    h2_order = order(h2s_results$posterior,decreasing = T)
    cum_posterior = cumsum(h2s_results$posterior[h2_order])
    
    num_top_grid = max(1,sum(cum_posterior <= target_prob))
    
    num_top_grid_marginal = sapply(1:(ncol(h2s_results)-1),function(i) {
      marginal_posterior = data.frame(h2 = h2s_results[,i],posterior = h2s_results$posterior)
      marginal_posterior = aggregate(posterior~h2,marginal_posterior,FUN = sum)
      max(1,sum(cumsum(sort(marginal_posterior$posterior,decreasing = T)) < target_prob))
    })
    
    top_h2s = h2s_results[h2_order[1:num_top_grid],1:ncol(tested_h2s),drop=FALSE]
    h2s_to_test = get_h2s_to_test(top_h2s,tested_h2s,step_size,ML,REML)
    rownames(h2s_to_test) = NULL
    
    if(nrow(h2s_to_test) == 0 && (num_top_grid < thresh_nonzero || min(num_top_grid_marginal) < thresh_nonzero_marginal)) {
      step_size = step_size/2
      h2s_to_test = get_h2s_to_test(top_h2s,tested_h2s,step_size,ML,REML)
      rownames(h2s_to_test) = NULL
    } 
    
    # if(sum(h2s_results$posterior[1:length(new_h2s_solutions)]) < 1-target_prob) break
  }
  h2s_results = h2s_results[do.call(order,lapply(seq_along(ncol(h2s_results)-1),function(x) h2s_results[,x])),,drop=FALSE]
  
  return(list(h2s_results = h2s_results,h2s_solutions = h2s_solutions,V_setup=V_setup))
}



GridLMM_ML = function(model,data,Y = NULL, X = NULL, weights = NULL,relmat = NULL,  
                     initial_step = 0.5,tolerance = 0.01,ML = T,REML=T,
                     V_setup = NULL, save_V_folder = NULL, # character vector giving folder name to save V_list
                     diagonalize=T,svd_K = T,drop0_tol = 1e-10,mc.cores = my_detectCores(),verbose=T) {
  
  # -------- check terms in models ---------- #
  terms = c(all.vars(model))
  if(!all(terms %in% colnames(data))) {
    missing_terms = terms[!terms %in% colnames(data)]
    stop(sprintf('terms %s missing from data',paste(missing_terms,sep=', ')))
  }
  
  # -------- Response ---------- #
  n = nrow(data)
  if(is.null(Y)){
    if(length(model) == 3){
      # if(length(model) == 2) then there is no response
      Y = as.matrix(data[,all.vars(model[[2]])])
      storage.mode(Y) = 'numeric'
      if(is.null(colnames(Y))) colnames(Y) = 1:ncol(Y)
    }
  }
  
  # -------- Constant Fixed effects ---------- #
  X_cov = model.matrix(nobars(model),data)
  if(!is.null(X)) X_cov = cbind(X_cov,X)
  linear_combos = caret::findLinearCombos(X_cov)
  if(!is.null(linear_combos$remove)) {
    cat(sprintf('dropping column(s) %s to make covariates full rank\n',paste(linear_combos$remove,sep=',')))
    X_cov = X_cov[,-linear_combos$remove]
  }
  if(any(is.na(X_cov))) stop('Missing values in covariates')
  
  # convert X_cov into X_list
  p = ncol(X_cov)
  
  # -------- Random effects ---------- #
  RE_setup = make_RE_setup(model = model,data = data,relmat = relmat)
  
  # -------- prep V ---------- #
  if(is.null(V_setup))  V_setup = make_V_setup(RE_setup,weights,diagonalize,svd_K,drop0_tol,save_V_folder,verbose)
  
  # -------- rotate data ---------- #
  # speeds to everything if a single random effect
  Qt = V_setup$Qt
  if(!is.null(Qt)){
    Y <- as.matrix(Qt %*% Y)
    X_cov <- as.matrix(Qt %*% X_cov)
  }
  
  # -------- prep_h2s ---------- #
  RE_names = names(RE_setup)
  n_RE = length(RE_names)
  
  step_size = initial_step
  h2s_to_test = make_h2s_matrix(RE_names,seq(0,1,by=step_size))
  
  # while the number of grid points accounting for target_prob % of the posterior is less than thresh_nonzero, keep dividing the h2s_matrix intervals smaller and smaller around the top regions
  h2s_solutions = c()
  tested_h2s = c()
  results = c()
  while(step_size >= tolerance && nrow(h2s_to_test)>0) {
    if(verbose) print(sprintf('step_size: %s, num_to_test: %d', step_size, nrow(h2s_to_test)))
    
    # -------- calculate likelihoods ----------- #
    registerDoParallel(mc.cores)
    results_list = foreach(h2s = iter(h2s_to_test,by = 'row')) %dopar% {
      chol_V_setup = make_chol_V_setup(V_setup,unlist(h2s))
      chol_V = chol_V_setup$chol_V
      V_log_det <- chol_V_setup$V_log_det
      inv_prior_X = rep(0,p)
      calc_LL(Y,X_cov,list(matrix(0,n,0)),t(h2s),chol_V,V_log_det,inv_prior_X,NULL,NULL,REML,BF=FALSE)
    } 
    tested_h2s = rbind(tested_h2s,h2s_to_test)
    if(length(results) > 0) {
      results = compile_results(c(list(results),results_list))
    } else {
      results = compile_results(results_list)
      step_size = step_size / 2
    }
    
    current_h2s = get_current_h2s(results,RE_names,ML,REML)
    h2s_to_test = get_h2s_to_test(current_h2s,tested_h2s,step_size,ML,REML)
    if(nrow(h2s_to_test) == 0 & step_size >= tolerance){
      step_size = step_size / 2
      h2s_to_test = get_h2s_to_test(current_h2s,tested_h2s,step_size,ML,REML)
    }
  }
  return(list(results = results, setup = V_setup))
}