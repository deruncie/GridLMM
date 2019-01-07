#' Evaluate the posterior of a linear mixed model using the Grid-sampling algorithm
#' 
#' Performs approximate posterior inference of a linear mixed model by evaluating the posterior over a grid
#'    of values of the variance component proportions. The grid evaluation uses a heuristic to avoid traversing 
#'    the entire grid. This should work well as long as the posterior is unimodal.        
#'
#' @details       
#' Posterior inference involves an adaptive grid search. Generally, we start with a very coarse grid (with as few as 2-3 vertices per variance component)
#'    and then progressively increase the grid resolution focusing only on regions of high posterior probability. This is controlled
#'    by \code{h2_divisions}, \code{target_prob}, \code{thresh_nonzero}, and \code{thresh_nonzero_matrginal}. The sampling algorithm is as follows:
#'    \itemize{
#'    \item Start by evaluating the posterior at each vertex of a trial grid with resolution \eqn{m}
#'    \item Find the minimum number of vertices needed to sum to \code{target_prob} of the current (discrete) posterior. 
#'       Repeat for the marginal posteriors of each variance component#'    
#'    \item If these numbers are smaller than \code{thresh_nonzero} or \code{thresh_nonzero_matrginal}, respectively, form a new grid
#'       by increasing the grid resolution to \eqn{m/2}. Otherwise, STOP.
#'    \item Begin evaluating the posterior at the new grid only at those grid vertices that are adjacent (in any dimension) to any of the top
#'       grid vertices in the old grid.
#'    \item Re-evaluate the distribution of the posterior over the new grid. If any new vertices contribute to the top \code{target_prob} fraction of the 
#'       overall posterior, include these in the "top" set and return to step 4. 
#'       Note - the prior weights for the grid vertices must be updated each time the grid increases in resolution.
#'    \item Repeat steps 4-5 until no new grid vertices contribute to the "top" set.
#'    \item Repeat steps 2-6 until a STOP is reached at step 3.
#'    }
#'    
#' Note: Default parameters for priors give flat (improper) priors. These should be used with care, especially for calculations of Bayes Factors.
#'
#' @param formula A two-sided linear formula as used in \code{\link[lme4]{lmer}} describing the fixed-effects
#'    and random-effects of the model on the RHS and the response on the LHS. Note: correlated random-effects
#'    are not implemented, so using one or two vertical bars (\code{|}) or one is identical. At least one random effect is needed.
#' @param data A data frame containing the variables named in \code{formula}.
#' @param weights An optional vector of observation-specific weights.
#' @param relmat A list of matrices that are proportional to the (within) covariance structures of the group level effects. 
#'     The names of the matrices should correspond to the columns in \code{data} that are used as grouping factors. All levels
#'     of the grouping factor should appear as rownames of the corresponding matrix.
#' @param h2_divisions Starting number of divisions of the grid for each variance component.
#' @param h2_prior Function that takes two arguments: 1) A vector of \code{h2s} (ie variance component proportions), and 
#'     2) An integer giving the number of vertices in the full grid. The function should return a non-negative value giving the prior 
#'     weight to the grid cell corresponding to \code{h2s}.
#' @param a,b Shape and Rate parameters of the Gamma prior for the residual variance \eqn{\sigma^2}. Setting both to zero gives a limiting "default" prior.
#' @param inv_prior_X Vector of values for the prior precision of each of the fixed effects (including an intercept). Will be recycled if necessary.
#' @param target_prob See \strong{Details}.
#' @param thresh_nonzero  See \strong{Details}.
#' @param thresh_nonzero_marginal  See \strong{Details}.
#' @param V_setup Optional. A list produced by a GridLMM function containing the pre-processed V decompositions for each grid vertex, 
#'     or the information necessary to create this. Generally saved from a previous run of GridLMM on the same data.
#' @param save_V_folder Optional. A character vector giving a folder to save pre-processed V decomposition files for future / repeated use. 
#'     If null, V decompositions are stored in memory
#' @param diagonalize If TRUE and the model includes only a single random effect, the "GEMMA" trick will be used to diagonalize V. This is done
#'     by calculating the SVD of K, which can be slow for large samples.
#' @param svd_K If TRUE and \code{nrow(K) < nrow(Z)}, then the SVD is done on \eqn{K} instead of \eqn{ZKZ^T}
#' @param drop0_tol Values closer to zero than this will be set to zero when forming sparse matrices.
#' @param mc.cores Number of cores to use for parallel evaluations.
#' @param verbose Should progress be printed to the screen?
#'
#' @return A list with three elements:
#' \item{h2s_results}{A data frame with each row an evaluated grid vertex, with the first \eqn{l} columns giving the \eqn{h^2}'s, and the final column the 
#'     corresponding posterior mass}
#' \item{h2s_solutions}{A list with the parameters of the NIG distribution for each grid vertex}
#' \item{V_setup}{The \code{V_setup} object for this model. Can be re-passed to this function (or other GridLMM functions) to re-fit the model to the same data.}
#' @export
#'
#' @examples
GridLMM_posterior = function(formula,data,weights = NULL,relmat = NULL,  
                  h2_divisions = 10,
                  h2_prior = function(h2s,n) 1/n, a = 0, b = 0, inv_prior_X = 0,
                  target_prob = 0.99, # want this much probability to be distributed over at least thresh_nonzero grid squares 
                  thresh_nonzero = 10, thresh_nonzero_marginal = 0,
                  V_setup = NULL, save_V_folder = NULL, # character vector giving folder name to save V_list
                  diagonalize=T,mc.cores = my_detectCores(),verbose=T) {
 
  MM = prepMM(formula,data,weights,other_formulas = NULL,
              relmat,X=NULL,X_ID=NULL,proximal_markers=NULL,V_setup,diagonalize, svd_K = TRUE,drop0_tol = 1e-10,save_V_folder, verbose)
  lmod = MM$lmod
  RE_setup = MM$RE_setup
  V_setup = MM$V_setup
  
  Y = matrix(lmod$fr[,1])
  colnames(Y) = 'y'
  X_cov = lmod$X
  data = lmod$fr
  
  p = ncol(X_cov)
  if(is.null(inv_prior_X)) {
    inv_prior_X = rep(0,p)
  } else if(length(inv_prior_X) == 1) {
    inv_prior_X = rep(inv_prior_X,p)
  } else if(length(inv_prior_X) != p) {
    stop('wrong length of inv_prior_X')
  }
  
  
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
      chol_Vi = chol_V_setup$chol_V
      V_log_det = 2*sum(log(diag(chol_Vi)))
      SSs <- GridLMM_SS_matrix(Y,chol_Vi,X_cov,NULL,integer(),inv_prior_X)
      
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
    
  }
  h2s_results = h2s_results[do.call(order,lapply(seq_along(ncol(h2s_results)-1),function(x) h2s_results[,x])),,drop=FALSE]
  
  return(list(h2s_results = h2s_results,h2s_solutions = h2s_solutions,V_setup=V_setup))
}



#' Find the (Restricted) Maximim Likelihood solution to a linear mixed model by grid search
#' 
#' Uses the GridLMM approach to find the ML or REML solution to variance componenents of a linear mixed model.
#' 
#' @details Note: this is not a particularly efficient technique for a single model. 
#'    However, the calculations for each grid vertex can be re-used in later fits (ex. GWAS).
#'    This function is used internally by other GridLMM functions to find starting values.
#'    
#' The function uses a hill-climbing heuristic on the grid. It starts be evaluating a grid with a step size of \code{initial_step} \eqn{h^2} units. 
#'     It then halves the step size and evaluates all grid vertices adjacent to the grid vertex with the highest ML (and/or REML) value.
#'     If a new maximum is found, the grid vertices adjacent to the new maximum are evaluated. This is repeated until no new maximum is found.
#'     At this point, if the step size is smaller than \code{tolerance}, the algorithm stops. Otherwise, the step size is halved again and the algorithm repeats.
#'
#' @inheritParams GridLMM_posterior
#' @param initial_step See \strong{Details}.
#' @param tolerance See \strong{Details}.
#' @param ML If TRUE, ML solutions will be found.
#' @param REML If TRUE, REML solutions will be found.
#'
#' @return A list with two elements:
#' \item{results}{A data frame with each row a ML or REML solution, and columns giving additional statistics and parameter values.}
#' \item{V_setup}{The \code{V_setup} object for this model. Can be re-passed to this function (or other GridLMM functions) to re-fit the model to the same data.}
#' @export
#'
#' @examples
GridLMM_ML = function(formula,data,weights = NULL,relmat = NULL,  
                     initial_step = 0.5,tolerance = 0.01,ML = T,REML=T,
                     V_setup = NULL, save_V_folder = NULL, # character vector giving folder name to save V_list
                     diagonalize=T,mc.cores = my_detectCores(),verbose=T) {
 
  V_setup = prepMM(formula,data,weights,other_formulas = NULL,
              relmat,X=NULL,X_ID=NULL,proximal_markers=NULL,V_setup,diagonalize, svd_K = TRUE,drop0_tol = 1e-10,save_V_folder, verbose)
  # lmod = MM$lmod
  # RE_setup = MM$RE_setup
  # V_setup = MM$V_setup
  
  Y = matrix(V_setup$fr[,1])
  colnames(Y) = 'y'
  X_cov = V_setup$X
  data = V_setup$fr
  p = ncol(X_cov)
  
  # -------- rotate data ---------- #
  # speeds up everything if a single random effect
  Qt = V_setup$reTrms$Qt
  if(!is.null(Qt)){
    Y <- as.matrix(Qt %*% Y)
    X_cov <- as.matrix(Qt %*% X_cov)
  }
  
  # -------- prep_h2s ---------- #
  RE_names = V_setup$reTrms$h2_names
  n_RE = length(RE_names)
  
  step_size = initial_step
  h2s_to_test = make_h2s_matrix(RE_names,seq(-1,1,by=step_size))
  h2s_to_test = check_h2s_matrix(h2s_to_test,V_setup$reTrms$is_var)
  
  # while the number of grid points accounting for target_prob % of the posterior is less than thresh_nonzero, keep dividing the h2s_matrix intervals smaller and smaller around the top regions
  h2s_solutions = c()
  tested_h2s = c()
  results = c()
  
  while(step_size >= tolerance && nrow(h2s_to_test)>0) {
    if(verbose) print(sprintf('step_size: %s, num_to_test: %d', step_size, nrow(h2s_to_test)))
    
    # -------- calculate likelihoods ----------- #
    registerDoParallel(mc.cores)
    results_list = foreach(h2s = iter(h2s_to_test,by = 'row')) %dopar% {
      # print(h2s)
      chol_V_setup = make_chol_V_setup(V_setup,unlist(h2s))
      if(is.null(chol_V_setup)) return()
      chol_Vi = chol_V_setup$chol_V
      inv_prior_X = rep(0,p)
      calc_LL(Y,X_cov,X_list=NULL,t(h2s),chol_Vi,inv_prior_X,NULL,NULL,REML,BF=FALSE)
      # SSs <- GridLMM_SS_matrix(Y,chol_Vi,X_cov,list(matrix(0,n,0)),integer(),inv_prior_X)
      # chol_V = chol_V_setup$chol_V
      # V_log_det <- chol_V_setup$V_log_det
      # calc_LL(Y,X_cov,list(matrix(0,n,0)),t(h2s),chol_V,V_log_det,inv_prior_X,NULL,NULL,REML,BF=FALSE)
    } 
    results_list = results_list[lengths(results_list)>0]
    tested_h2s = rbind(tested_h2s,h2s_to_test)
    if(length(results) > 0) {
      results = compile_results(c(list(results),results_list))
    } else {
      results = compile_results(results_list)
      step_size = step_size / 2
    }
    
    current_h2s = get_current_h2s(results,RE_names,ML,REML)
    h2s_to_test = get_h2s_to_test(current_h2s,tested_h2s,step_size,ML,REML,V_setup$reTrms$is_var)
    if(nrow(h2s_to_test) == 0 & step_size >= tolerance){
      step_size = step_size / 2
      h2s_to_test = get_h2s_to_test(current_h2s,tested_h2s,step_size,ML,REML,V_setup$reTrms$is_var)
    }
    # recover()
  }
  return(list(results = results, setup = V_setup))
}