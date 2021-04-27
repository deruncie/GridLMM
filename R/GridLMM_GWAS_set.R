
#' GridLMM GWAS set
#'
#' Performs a GWAS of set-tests using GridLMM algorithm. Can perform LRTs, Wald-tests, or calculate Bayes Factors.
#' By default, uses the targeted grid search heuristic (fast algorithm), though can perform a full grid search as well.
#'
#' @details GridLMM performs approximate likelihood or posterior-based inference for linear mixed models efficiently by
#'  finding solutions to many models in parallel. Rather than optimizing to high precision for each separate model, GridLMM
#'  finds "good enough" solutions that satisfy many tests at once - so the expensive calculations can be re-used. It does
#'  this by trying variance components on a grid, and selecting the best grid cell for each model. The \code{Full} algorithm
#'  performs a full grid search over all variance component parameters. The \code{Fast} algorithm uses heuristics to reduce
#'  the number of grid cells that need to be evaluated - focusing from the maximum likelihood solutions under a null model
#'  with no markers, and then working out to neighboring grid cells from there.
#'
#'  Posterior inference involves an adaptive grid search. Generally, we start with a very coarse grid (with as few as 2-3 vertices per variance component)
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
#' @param formula A two-sided linear formula as used in \code{\link[lme4]{lmer}} describing the fixed-effects
#'    and random-effects of the model on the RHS and the response on the LHS. Note: correlated random-effects
#'    are not implemented, so using one or two vertical bars (\code{|}) or one is identical. At least one random effect is needed.
#'    Unlike \code{lmer}, random effects can have as many as there are observations.
#' @param test_formula test_formula One-sided formula for the alternative model (ML or BF), or full model (REML) to be applied to each test (ie \emph{marker}, or column of \code{X}).
#'     Each term on the RHS will be multiplied by a column of X to form a new covariate.
#'     Ex. \code{~1} specifices an intercept for each marker.
#'     \code{~1+cov} species an intercept and slope on \code{cov} for each marker.
#' @param reduced_formula One-sided formula for the reduced model. Same format as \code{test_formula}. Should have fewer degrees of freedom than \code{test_formula}. Not used for REML models.
#' @param data A data frame containing the variables named in \code{formula}.
#' @param weights An optional vector of observation-specific weights.
#' @param X Matrix of markers with \eqn{p} columns.
#'     Each column of X is used as a separate association test.
#'     Should have row names that correspond to the \code{X_ID} column of \code{data}.
#'     Colnames are used as IDs for each test, and should align with names of \code{proximal_markers} if provided.
#' @param X_ID Column of \code{data} that identifies the row of \code{X} that corresponding to each observation.
#'     It is possible that multiple observations reference the same row of \code{X}.
#' @param centerX,scaleX,fillNAX TRUE/FALSE for each. Applied to the \code{X} matrix before using \code{X} to form any GRMs.
#' @param X_map Optional. Data frame with information on each marker such as chromosome, position, etc. Will be appended to the results
#' @param relmat Either:
#'     1) A list of matrices that are proportional to the (within) covariance structures of the group level effects.
#'     2) A list of lists with elements (\code{K}, \code{p}) with a covariance matrix and an integer listing the number of markers
#'         used to estimate the covariance matrix. This is used for appropriate downdating of \code{V} to remove proximal markers for each test.
#'     The names of the matrices / list elements should correspond to the columns in \code{data} that are used as grouping factors. All levels
#'     of the grouping factor should appear as rownames of the corresponding matrix.
#' @param h2_step Step size of the grid
#' @param h2_start Optional. Matrix with each row a vector of \eqn{h^2} parameters defining starting values for the grid.
#'     Typically ML/REML solutions for the null model.
#'     If null, will be calculated using \link{GridLMM_ML}.
#' @param h2_start_tolerance Optional. Grid size for \link{GridLMM_ML} in finding ML/REML solutions for the mull model.
#' @param max_steps Maximum iterations of the heuristic algorithm per marker.
#' @param method One of 'REML', 'ML', or 'BF'.
#'     'REML' wimplies a Wald-test. 'ML' implies Maximum Likelihood evaluation, with the LRT.
#'     'BF' does posterior evaluation and calculates Bayes Factors.
#' @param algorithm Either 'Fast' or 'Full'. See details.
#' @param inv_prior_X Vector of values for the prior precision of each of the fixed effects (including an intercept). Will be recycled if necessary.
#' @param target_prob see \strong{Details}
#' @param proximal_markers A list of integer vectors with length equal to the number of columns of \code{X}/
#'     Each element is a vector of indices of markers that should be removed from any GRMs before the current test is calculated.
#'     If \code{proximal_Xs} is provided, then the indices correspond to columns of \code{proximal_Xs}.
#'     Otherwise, the indices correspond to columns of \code{X}.
#'     If null, no downdating will be performed.
#' @param proximal_Xs Optional. A list of matrices to be used for downdating GRMs. If multiple GRMs are calculated from markers,
#'     this list can have multiple elements. Each matrix should have rownames like \code{X} corresponding to the levels of \code{X_ID} in \code{data}.
#' @param V_setup Optional. A list produced by a GridLMM function containing the pre-processed V decompositions for each grid vertex,
#'     or the information necessary to create this. Generally saved from a previous run of GridLMM on the same data.
#' @param save_V_folder Optional. A character vector giving a folder to save pre-processed V decomposition files for future / repeated use.
#'     If null, V decompositions are stored in memory
#' @param diagonalize If TRUE and the model includes only a single random effect, the "GEMMA" trick will be used to diagonalize V. This is done
#'     by calculating the SVD of K, which can be slow for large samples.
#' @param mc.cores Number of processor cores used for parallel evaluations. Note that this uses 'mclapply', so the memory requires grow rapidly with \code{mc.cores}, because
#'     the marker matrix gets duplicated in memory for each core.
#' @param verbose Should progress be printed to the screen?
#'
#'
#' @return A list with two elements:
#' \item{results}{A data frame with each row the results of the association test for a column of \code{X}, plus asssociated parameter values and statistics.}
#' \item{setup}{A list with several objects needed for re-running the model, including \code{V_setup} and \code{downdate_Xs}.
#'   These can be re-passed to this function (or other GridLMM functions) to re-fit the model to the same data.}
#' @export
#'
GridLMM_GWAS_set = function(formula,data,weights = NULL,
                        X,X_ID = 'ID', set_matrix,
                        centerX = FALSE,scaleX = FALSE,fillNAX = FALSE,X_map = NULL, relmat = NULL, normalize_relmat = TRUE,
                        h2_step = 0.01, h2_start = NULL, h2_start_tolerance = 0.001,max_steps = 100,
                        method = c('REML'), algorithm = c('Fast','Full'),
                        inv_prior_X = NULL,target_prob = 0.99,
                        proximal_markers = NULL, proximal_Xs = NULL,
                        V_setup = NULL, save_V_folder = NULL,
                        diagonalize=T,
                        mc.cores = my_detectCores(),
                        verbose=T) {

  # Evaluate argument options
  # -------- Evaluate argument options ---------- #
  method = match.arg(method)
  algorithm = match.arg(algorithm)

  
  # -------- check set_matrix ---------- #
  if(!is.matrix(set_matrix)) set_matrix = as.matrix(set_matrix)
  if(ncol(set_matrix) != 2) stop("set_matrix must have 2 columns")
  if(min(set_matrix) < 1 || max(set_matrix) > ncol(X)) stop("set_matrix refers to columns not in X")
  if(is.null(rownames(set_matrix))) rownames(set_matrix) = paste0('Set_',1:nrow(set_matrix))
  
  # -------- Re-scale SNPs ---------- #
  X = scale_SNPs(X,centerX,scaleX,fillNAX)

  # -------- prep Mixed Models ---------- #
  MM = prepMM(formula,data,weights,other_formulas = NULL,
              relmat,normalize_relmat,X,X_ID,proximal_markers,V_setup,diagonalize, svd_K = TRUE,drop0_tol = 1e-10,save_V_folder, verbose)
  lmod = MM$lmod
  RE_setup = MM$RE_setup
  V_setup = MM$V_setup

  Y = matrix(lmod$fr[,1])
  colnames(Y) = 'y'
  X_cov = lmod$X
  data = lmod$fr

  if(is.null(colnames(X))) colnames(X) = paste('X',1:ncol(X),sep='_')

  # -------- prep Markers and proximal_Xs ---------- #
  # check that proximal_Xs is appropriate
  if(!is.null(proximal_markers)) {
    if(length(proximal_markers) != nrow(set_matrix)) stop("If provided, proximal_markers must have length == nrow(set_matrix)")
    if(is.null(names(proximal_markers))) names(proximal_markers) = rownames(set_matrix)
    if(length(RE_setup) > 1) {
      # downdate_REs is a vector of length == length(RE_setup) with 1s identifying REs that will be downdated
      downdate_REs = sapply(names(RE_setup),function(x) strsplit(x,'.',fixed=TRUE)[[1]]) == X_ID
      downdate_REs = as.numeric(downdate_REs)

      if(is.null(proximal_Xs)) {
        # use X_list_full
        if(length(downdate_REs) != 1) stop("Can't tell what markers to downdate from each RE")
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
# 
#   # -------- set up prior for X ------ #
#   if(method == 'BF') {
#     if(is.null(inv_prior_X)) {
#       warning('Calculating Bayes Factors with impropper prior on X might be unreliable')
#       inv_prior_X = rep(0,ncol(X_cov) + ncol(mm_test))
#     } else if(length(inv_prior_X) == 1) {
#       inv_prior_X = c(rep(0,ncol(X_cov)),rep(inv_prior_X,ncol(mm_test)))
#     } else if(length(inv_prior_X) == ncol(mm_test)) {
#       inv_prior_X = c(rep(0,ncol(X_cov)),inv_prior_X)
#     } else if(length(inv_prior_X) != (ncol(X_cov)+ncol(mm_test))) {
#       stop('Wrong length of inv_prior_X')
#     }
#   } else{
#     inv_prior_X = rep(0,ncol(X_cov) + ncol(mm_test))
#   }

  # -------- get starting value under null model ---------- #
  # evaluate likelihoods on a grid, going until all LL's stop increasing
  if(!is.null(Y) && is.null(h2_start)) {# && method != 'BF') {
    if(verbose) print('Estimating h2_start via null model')
    if(ncol(Y) > 1) stop('null model h2 only implemented for a single response')
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
    h2_start = t(setup_Grid(names(RE_setup),h2_step,h2_start))
  }


  setup = list(
    Y = Y,
    X_cov = X_cov,
    set_matrix = set_matrix,
    h2_start = h2_start,
    V_setup = V_setup,
    proximal_markers = proximal_markers,
    proximal_Xs = proximal_Xs
  )


  results = run_GridLMM_GWAS_set(Y,X_cov,X,set_matrix,inv_prior_X,X_map,V_setup,h2_start,h2_step,target_prob,proximal_markers, proximal_Xs, method,mc.cores,verbose)

  return(list(
    results = results,
    setup = setup
  ))
}


run_GridLMM_GWAS_set = function(Y,X_cov,X,set_matrix,inv_prior_X = NULL,X_map = NULL, V_setup,h2_start=NULL,h2_step=0.01,target_prob = 0.99,
                            proximal_markers=NULL, proximal_Xs=NULL,
                            method = c('REML'), mc.cores=my_detectCores(),verbose = FALSE)
{

  # -------- Evaluate argument options ---------- #
  method = match.arg(method)

  # -------- check inputs ------ #
  n = nrow(Y)
  if(nrow(X_cov) != n) stop("Wrong dimensions of X_cov")
  if(nrow(X) != n) stop("Wrong dimensions of X")
  if(is.null(rownames(set_matrix))[1]) rownames(set_matrix) = paste0("Set_",1:nrow(set_matrix))
  if(!is.null(proximal_markers) && is.null(names(proximal_markers)[1])) names(proximal_markers) = rownames(set_matrix)

  if(is.null(inv_prior_X)) {
    if(method == 'BF') stop("Can't calculate Bayes Factors with impropper prior on X")
    inv_prior_X = rep(0,ncol(X_cov))
  }
  if(!length(inv_prior_X) == (ncol(X_cov))) stop("Wrong length of inv_prior_X")


  # -------- rotate data ---------- #
  # speeds to everything if a single random effect

  Qt = V_setup$Qt
  if(!is.null(Qt)){
    Y <- as.matrix(Qt %*% Y)
    X_cov <- as.matrix(Qt %*% X_cov)
    X <- as.matrix(Qt %*% X)
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

  results = fit_GridLMM_GWAS_set(Y,X_cov,X,set_matrix,h2_start,h2_step,V_setup,inv_prior_X,target_prob,proximal_markers,proximal_Xs,method,verbose,mc.cores)
  
  # Do tests
  # if(method == 'REML') {
  #   F_hat_cols = colnames(results)[grep('F.',colnames(results),fixed=T)]
  #   p_value_REML = do.call(cbind,lapply(F_hat_cols,function(col) p_value = pf(results[,col],1,nrow(Y) - results$Df_X - ncol(X_cov),lower.tail=F)))
  #   if(length(F_hat_cols) > 1) {
  #     colnames(p_value_REML) = sprintf('p_value_REML.%d',1:length(F_hat_cols))
  #   } else{
  #     colnames(p_value_REML) = 'p_value_REML'
  #   }
  #   results = data.frame(results,p_value_REML)
  # } else {
    # Need to provide proximal_Xs for reduced model if proximal_Xs draws from X_list_full
  if(is.null(proximal_Xs) || !is.list(proximal_Xs)) {
    if(!is.list(proximal_Xs)) {
      proximal_Xs_new = list()
      proximal_Xs_new[[proximal_Xs]] = X
      proximal_Xs = proximal_Xs_new
      rm(X)
      gc()
    }
  }
  if(method == 'REML') {
    results_reduced = fit_GridLMM_GWAS_set(Y,X_cov,NULL,set_matrix,h2_start,h2_step,V_setup,inv_prior_X,target_prob,proximal_markers,proximal_Xs,method,verbose)
    for(trait in unique(results$Trait)){
      index_full = results$Trait == trait
      index_reduced = results_reduced$Trait == trait
      results$REML_Reduced_logLik[index_full] = results_reduced$REML_logLik[index_reduced]
      # results$Reduced_Df_X[index_full] = results_reduced$Df_X[index_reduced]
      # results$Reduced_X_id[index_full] = results_reduced$X_id[index_reduced]

      # results$p_value_ML = with(results,pchisq(2*(ML_logLik - ML_Reduced_logLik),Df_X - Reduced_Df_X,lower.tail = F))
    }
  }
    # if(method == 'BF'){
    #   results_reduced = fit_GridLMM_GWAS(Y,X_cov,X_list_reduced,h2_start,h2_step,V_setup,inv_prior_X,target_prob,proximal_markers,proximal_Xs,method,verbose)
    #   for(trait in unique(results$Trait)){
    #     index_full = results$Trait == trait
    #     index_reduced = results_reduced$Trait == trait
    #     results$log_posterior_factor_reduced = results_reduced$log_posterior_factor[index_reduced]
    #     results$BF[index_full] = results$log_posterior_factor - results$log_posterior_factor_reduced # Note: need to account for different V_beta
    #     results$Reduced_Df_X[index_full] = results_reduced$Df_X[index_reduced]
    #     if(!is.null(inv_prior_X)){
    #       if(length(inv_prior_X) == 1) {
    #         results$BF[index_full] = results$BF[index_full] - log(inv_prior_X)/2
    #       } else{
    #         results$BF[index_full] = results$BF[index_full] - sum(log(tail(inv_prior_X,n = results$Df_X[1] - results$Reduced_Df_X[1])))/2
    #       }
    #     }
    #   }
    # }
  # }
  results$ML_index = c()
  results$REML_index = c()

  # add map info
  if(!is.null(X_map)){
    index = match(results$X_ID,X_map$snp)
    if(!is.null(index) && length(index) == nrow(results)) {
      results = data.frame(X_map[index,],results)
    }
  }

  return(results)
}


fit_GridLMM_GWAS_set = function(Y,X_cov,X,set_matrix,h2_start,h2_step,V_setup,inv_prior_X = NULL,target_prob = 0.99,
                            proximal_markers = NULL,proximal_Xs = NULL,
                            method = 'ML',verbose=F,mc.cores = 1) {
  # calculates LL's for each test over a grid of h2s
  # starts at h2_start, and then moves out in equidistant circles.
  # Only tests where the LL increases in one circle are re-tested in the next.
  # continues until the grid is fully covered, or no tests increase in LL
  Y = as.matrix(Y)
  n = nrow(Y)
  storage.mode(Y) = 'double'
  if(is.null(colnames(Y))) colnames(Y) = 1:ncol(Y)
  
  if(ncol(set_matrix) == 2) {
    set_matrix = cbind(1:nrow(set_matrix),set_matrix)
  } else{
    stop("set_matrix must have 2 columns")
  }

  if(method != 'REML') stop("Only REML implemented")
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
  if(is.null(X)) {
    X = matrix(0,n,0)
    if(is.null(proximal_markers)) {
      null_model = TRUE
      stop("model without proximal_markers not implemented. Can be an empty list of length == nrow(set_matrix)")
    } else{
      if(is.null(names(proximal_markers))) stop("proximal_markers must have names")
    }
  } else if(is.null(colnames(X)) & ncol(X)>0) colnames(X) = 1:ncol(X)

  if(is.null(inv_prior_X)){
    inv_prior_X = rep(0,ncol(X_cov))
  }

  RE_names = colnames(h2_start)
  if(is.null(RE_names) && !is.null(names(h2_start))) RE_names = names(h2_start)
  if(null_model) {
    # active = TRUE
    # h2s_to_test = as.matrix(h2_start)
    # tested_h2s = c()
    # results = c()
    # 
    # n_steps = 1
    # while(active && nrow(h2s_to_test) > 0) {
    #   if(verbose) print(sprintf('step: %d, num_h2s: %d, num active: %d',n_steps,nrow(h2s_to_test), sum(active)))
    #   inner_cores = 1
    #   outer_cores = mc.cores
    #   registerDoParallel(outer_cores)
    #   results_list = foreach(h2s = iter(h2s_to_test,by='row')) %dopar% {
    #     h2s = h2s[1,]
    #     chol_V_setup = make_chol_V_setup(V_setup,h2s)
    #     chol_V = chol_V_setup$chol_V
    #     calc_LL_parallel(Y,X_cov,X_list,h2s,chol_V,inv_prior_X[1:ncol(X_cov)],
    #                      NULL,V_setup$n_SNPs_downdated_RRM,REML,BF,inner_cores,c())
    #   }
    #   tested_h2s = rbind(tested_h2s,h2s_to_test)
    #   if(length(results) == 0) {
    #     results = compile_results(results_list)
    #   } else{
    #     results = compile_results(c(list(results),results_list))
    #   }
    # 
    #   if(n_steps > 1) {
    #     active = FALSE
    #     if(ML && !is.null(results$ML_index) && results$ML_index > 1) active = TRUE
    #     if(REML && !is.null(results$REML_index) && results$REML_index > 1) active = TRUE
    #     if(BF && results$BF_new_posterior_mass > 1-target_prob) active = TRUE
    #   }
    #   if(active) {
    #     results$n_steps = n_steps
    #   }
    # 
    #   if(!BF) {
    #     current_h2s = get_current_h2s(results,RE_names,ML,REML)
    #   } else{
    #     current_h2s = tested_h2s
    #   }
    #   h2s_to_test = get_h2s_to_test(current_h2s,tested_h2s,h2_step,ML,REML)
    #   if(ncol(h2s_to_test) == 0) break
    #   n_steps = n_steps + 1
    # }
    # 
    # results$Df_X = rep(0,ncol(Y))

  } else{

    h2s_to_test = as.matrix(h2_start)
    tested_h2s = c()

    if(is.null(proximal_markers)) {
      active_tests = rep(TRUE,ncol(X))
    } else{
      active_tests = rep(TRUE,length(proximal_markers))
    }
    n_steps = 1
    results = c()

    while(sum(active_tests) > 0 && nrow(h2s_to_test) > 0) {
      if(verbose) print(sprintf('step: %d, num_h2s: %d, num active: %d',n_steps,nrow(h2s_to_test), sum(active_tests)))
      active_set_matrix = set_matrix[active_tests,,drop=FALSE]

      if(!is.null(proximal_markers)) {
        if(is.list(proximal_Xs)){
          downdate_Xs = build_downdate_Xs(1:length(proximal_Xs),proximal_Xs,proximal_markers,active_set_matrix[,1])
        } else{
          downdate_Xs = build_downdate_Xs(proximal_Xs,list(X),proximal_markers,active_set_matrix[,1])
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
        calc_LL_set(Y,X_cov,X,h2s,chol_V,inv_prior_X,
                         downdate_Xs,V_setup$n_SNPs_downdated_RRM,REML,BF,active_set_matrix)
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
        # if(ML && !is.null(results$ML_index)) active_tests[!is.na(results$ML_index) & results$ML_index > 1] = TRUE
        if(REML && !is.null(results$REML_index)) active_tests[!is.na(results$REML_index) & results$REML_index > 1] = TRUE
        # if(BF) active_tests[results$BF_new_posterior_mass > 1-target_prob] = TRUE
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
    # results$Df_X = rep(length(X_list),ncol(Y))
    # if(length(X_list) == 0) results$Df_X = 0

  }
  return(results)
}


calc_LL_set = function(Y,X_cov,X,h2s,chol_Vi,inv_prior_X,
                   downdate_Xs = NULL,n_SNPs_downdated_RRM = NULL,REML = TRUE, BF = TRUE,active_set_matrix){  
  
  m = ncol(Y)
  n = nrow(Y)
  if(is.null(downdate_Xs)) {
    stop("no downdate_Xs not implemented")
    # SSs <- GridLMM_SS_matrix(Y,chol_Vi,X_cov,X_list,active_X_list,inv_prior_X)
  } else{
    downdate_weights = h2s/unlist(n_SNPs_downdated_RRM)
    if(length(downdate_weights) != length(downdate_Xs[[active_set_matrix[1,1]]]$downdate_Xi)) stop("Wrong length of downdate weights")
    chol_Vi = as.matrix(chol_Vi)
    SSs <- GridLMM_setTest_downdate_matrix(Y,chol_Vi,X_cov,X,active_set_matrix,downdate_Xs,downdate_weights,inv_prior_X);
  }
  log_LL = SSs
  # log_LL = get_LL(SSs,X_cov,X_list,active_X_list,n,m,ML=TRUE,REML = REML,BF=BF)
  
  # collect results in a table
  # beta_hats and F_hats are lists of matrices, with coefficients as columns and traits as rows. 
  # We want to organize into a single matrix with coefficients as columns and each trait sequentially, with all tests for each trait consecutive
  
  # p = max(1,max(ncol(X_list[[1]]),length(downdate_weights_i)))*m
  # ML_h2 = as.matrix(t(h2s))[rep(1,p),,drop=FALSE]
  n_tests = length(log_LL$REML)
  if(is.null(n_tests)) n_tests = 1
  p = n_tests*m
  # ML_h2 = as.matrix(t(h2s))[rep(1,p),,drop=FALSE]
  # colnames(ML_h2) = paste(colnames(ML_h2),'ML',sep='.')
  
  results_i = data.frame(Trait = rep(colnames(Y),each = n_tests),
                         X_ID = NA,
                         # ML_logLik = c(log_LL$ML),
                         # ML_h2,
                         stringsAsFactors = F)
  if(nrow(active_set_matrix) == n_tests) {
    if(ncol(X)>0) {
      results_i$X_ID = colnames(X)[active_set_matrix[,1]]
    } else {
      results_i$X_ID = names(downdate_Xs)[active_set_matrix[,1]]
    }
  } else{
    results_i$X_ID = 'NULL'
  }
  
  b_cov = ncol(X_cov)
  # b_x = length(X_list)
  # if(length(X_list) == 0) b_x = 0
  b = b_cov
  # if(ncol(X_list[[1]]) == n_tests) {
  # beta_hats = do.call(rbind,lapply(1:m,function(i) t(log_LL$beta_hat[(i-1)*b + b_cov + 1:b_x,,drop=FALSE])))
  beta_hats = do.call(rbind,lapply(1:m,function(i) t(log_LL$beta_hat[(i-1)*b + 1:b,,drop=FALSE])))
  colnames(beta_hats) = paste('beta',1:ncol(beta_hats),sep='.')
  results_i = data.frame(results_i,beta_hats)
  # }
  
  # if(REML) {
    REML_h2 = as.matrix(t(h2s))[rep(1,p),,drop=FALSE]
    colnames(REML_h2) = paste(colnames(REML_h2),'REML',sep='.')
    results_i = data.frame(results_i,REML_logLik = c(log_LL$REML),REML_h2,Tau = log_LL$tau)
    # if(b_x > 0) {
    #   F_hats = do.call(rbind,lapply(1:m,function(i) t(log_LL$F_hat[(i-1)*b + b_cov + 1:b_x,,drop=FALSE])))
    #   colnames(F_hats) = paste('F',1:ncol(F_hats),sep='.')
    #   results_i = data.frame(results_i,F_hats)
    # }
  # }
  # if(BF){
  #   results_i$log_posterior_factor = c(t(log_LL$log_posterior_factor))
  #   results_i$RSSs = c(t(log_LL$RSSs))
  # }
  return(results_i)
}

