#' LASSO solutions in a linear mixed model using GridLMM
#'
#' Finds LASSO or Elastic Net solutions for a multiple regression problem with correlated errors.
#' 
#' @details Finds the full LASSO or Elastic Net solution path by running \code{\link[glmnet]{glmnet}} at each grid vertex
#'
#' @inheritParams GridLMM_GWAS
#' @inheritParams glmnet::glmnet
#' @param X Variables in model that well be penalized with the elastic net penalty. Covariates specified in \code{formula} are not penalized.
#' @param ... 
#'
#' @return An object with S3 class "glmnet","*" , where "*" is "elnet". See \code{\link[glmnet]{glmnet}}.
#' @export
#'
#' @examples
GridLMMnet = function(formula,data,X, weights = NULL, 
                     centerX = TRUE,scaleX = TRUE,relmat = NULL,
                     h2_divisions = 10, h2_start = NULL,
                     alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL,
                     RE_setup = NULL, V_list_setup = NULL, save_V_list = NULL,
                     diagonalize=T,svd_K = T,drop0_tol = 1e-10,mc.cores = parallel::detectCores(),clusterType = 'mclapply',verbose=T,...) 
{
  
  # -------- check terms in formulas ---------- #
  terms = c(all.vars(formula))
  if(!all(terms %in% colnames(data))) {
    missing_terms = terms[!terms %in% colnames(data)]
    stop(sprintf('terms %s missing from data',paste(missing_terms,sep=', ')))
  }
  
  # -------- Response ---------- #
  n = nrow(data)
  if(length(formula) == 3){
    y = as.matrix(data[,all.vars(formula[[2]])])
  } else{
    stop('Check if formula missing LHS or RHS.')
  }
  nobs = n
  
  # -------- Constant Fixed effects ---------- #
  X_cov = model.matrix(nobars(formula),data)
  linear_combos = caret::findLinearCombos(X_cov)
  if(!is.null(linear_combos$remove)) {
    cat(sprintf('dropping column(s) %s to make covariates full rank\n',paste(linear_combos$remove,sep=',')))
    X_cov = X_cov[,-linear_combos$remove]
  }
  if(any(is.na(X_cov))) stop('Missing values in covariates')
  
  nvars = ncol(X_cov) + ncol(X)
  
  # -------- Random effects ---------- #
  if(is.null(RE_setup)) {
    RE_setup = make_RE_setup(formula = formula,data,relmat = relmat,verbose=verbose)
  }
  
  V_list_setup = make_V_list(RE_setup,
                             weights,
                             h2_divisions,  
                             h2_start,
                             save_V_list,
                             diagonalize,
                             svd_K,
                             drop0_tol,
                             mc.cores,
                             clusterType,
                             verbose)
  
  setup = list(
    X_cov = X_cov,
    alpha = alpha,
    V_list_setup = V_list_setup
  )
  if(is.null(y)) return(setup)
  
  results = GridLMMnet2(y,X_cov,X,V_list_setup,mc.cores,clusterType=clusterType,alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda=lambda,verbose=verbose,...)
  
  # returning elnet object. Don't think it's worth including the setup object too. I can add functionality to re-use V_list if passed.
  return(results)
  
  # return(list(
  #   res.glmnet = results,
  #   setup = setup
  # ))
}

run_glmnet_V = function(X_full,y,chol_V_setup, alpha = alpha, 
                      lambda = NULL,nlambda, lambda.min.ratio,get_max_lambda = FALSE,penalty.factor,...){
  nobs = length(y)
  
  if(is.character(chol_V_setup)) {
    chol_V_setup <- readRDS(chol_V_setup)
  }
  h2_i <- chol_V_setup$h2_index
  chol_V = chol_V_setup$chol_V
  V_log_det <- chol_V_setup$V_log_det
  
  chol_V_inv = solve(chol_V)
  
  # calculate a normalization constant based on V_inv
  # ci = exp(-V_log_det/nobs)  # normalization constant based on determinant
  ci = sum(chol_V_inv^2)/nobs # normalization constant based on trace of V_inv
  
  # normalize chol_V_inv
  chol_V_inv = chol_V_inv/sqrt(ci)
  V_log_det = -2*sum(log(diag(chol_V_inv)))
  
  y_star = as.matrix(crossprod(chol_V_inv,y))
  X_star = as.matrix(crossprod(chol_V_inv,X_full))
  sd_y_star = sd(y_star)*(nobs-1)/nobs
  
  intercept = FALSE
  if(var(X_star[,1]) == 0) intercept = TRUE
  
  if(get_max_lambda) {
    # code from: https://stats.stackexchange.com/questions/166630/glmnet-compute-maximal-lambda-value?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    # to calculate the maximum lambda without calling glmnet
    # modified to exclude intercept
    #    this is probably not quite right with non-penalized coefficients.
    #    probably need to actually project out these coefficients from X_star, or something like that.
    #    may be ok as an approximation
    #    also, can just run glmnet and take the value
    
    # Notes 6/26/18
    #    penalty factors are internally rescaled to sum to nvars, and the lambda sequence will reflect this change.
    #    max(penalty.factor*sum(penalty.factor)/(length(penalty.factor))*abs(t(resid(lm(y~X[,penalty.factor==0]))) %*% X ) / ( alpha * n))
    #    the above code seems to work, and is consistent with the paper
    #    Need to re-derive if/how lambda should change depending on sd(y_star)
    #    I wonder if modifiying penalty.factor as above will make the deviance calculation consistent?
    resid_y_star = resid(lm(y_star~X_star[,penalty.factor == 0]))
    s2_hat = sum(resid_y_star^2)/nobs
    max_lambda = max(penalty.factor*sum(penalty.factor)/(length(penalty.factor)-1)*abs(t(resid_y_star) %*% X_star ) / ( alpha * nobs))
    score = nobs*log(2*pi) + nobs*log(s2_hat) + V_log_det + nobs
    return(c(max_lambda=max_lambda,score=score))
  }
  
  # notes:
  # intercept = FALSE is needed, because intercepts are included in X_cov. However, if all(X_cov[,1] == 1), glmnet seems to ignore this.
  #      (Actually, it ignores any covariates with variance = 0)
  #   in this case, use the original intercept, but in results, put coefficient back into place. This requires that the intercept is the first column of X_cov
  # standardize = FALSE. Necessary because X_full is modified by X_star. We don't want it standardized differently for each Vi
  # standardize.response = FALSE. Want this, but only works with family='mgaussian'. 
  #   Instead, we first calculate a lambda sequence that should work for all Vi
  #   Then, when we use the lambda sequence, we normalize lambda, y, and X all by sd(y). This means normalizing lambda by sd(y)^2
  #   This way, the actual magnitude of y_star matters, and the results are comparable. 
  #   And the total scores can be calculate correctly.
  
  if(is.null(lambda)) {
    # run with standardized data, then adjust lambda for V_1.
    res = glmnet(X_star/sd_y_star,y_star/sd_y_star,family = 'gaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,
                 alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,intercept = intercept,standardize = FALSE,...)
    lambda = res$lambda * sd_y_star^2
    # alternatively, use family = 'mgaussian' and standarize.response = F
    # res = glmnet(X_star,y_star,family = 'mgaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,standardize.response = F,
    #              alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,intercept = intercept,standardize = FALSE,...)
    # class(res)[1] = 'elnet' # so we can use plot.elnet
    # lambda = res$lambda
  } else{
    # lambda = lambda * ci
    res = glmnet(X_star/sd_y_star,y_star/sd_y_star,family = 'gaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,
                 alpha = alpha, lambda = lambda/sd_y_star^2 ,intercept = intercept,standardize = FALSE,...)
    # alternatively, use family = 'mgaussian' and standarize.response = F
    # res = glmnet(X_star,y_star,family = 'mgaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,standardize.response = F,
    #               alpha = alpha, lambda = lambda,intercept = intercept,standardize = FALSE,...)
    # class(res)[1] = 'elnet' # so we can use plot.elnet
  }
  if(intercept){
    res$beta[1,] = res$a0*sd_y_star
    res$a0[] = 0
  }
  RSSs = colSums((y_star[,1] - as.matrix(X_star %*% res$beta))^2)
  penalty.factor = penalty.factor*(length(penalty.factor)-1)/sum(penalty.factor) # adjust penalty.factor to sum to p-1 (for intercept)
  penalties = colSums(sweep(penalty.factor*((1-alpha)/2*res$beta^2 + alpha*abs(res$beta)),2,lambda,'*')) # from: https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
  
  # score function taken from: Schelldorfer, J., Bühlmann, P., & DE GEER, S. V. (2011). Estimation for High-Dimensional Linear Mixed-Effects Models Using ℓ1-Penalization. Scandinavian Journal of Statistics, 38(2), 197–214. http://doi.org/10.1111/j.1467-9469.2011.00740.x
  # except that their V includes \hat{sigma}. This doesn't make sense to me. This parameter isn't in the normal LASSO score.
  # I've taken it out here, and it seems to work, I think.
  # Rakitsch et al 2013 Bioinformatics also shows the LASSOlmm score, but they have a fixed residual covariance, and so doesn't have to worry about |V|.
  # scores = V_log_det/2 + RSSs/2 + nobs*penalties  # multiply penalties by nobs because of how glmnet calculates penalty. #not sure if this is needed: nobs*log(s2_hats)/2 + 
  
  # -2 times penalized log-Likelihood function. This gives basically the same different answers. Slightly different h2s.
  # s2_hats = RSSs/nobs
  # scores = (nobs*log(2*pi) + nobs*log(s2_hats) + V_log_det + nobs) + 2*nobs*penalties
  
  # I think this one is correct. I've changed the log-likelihood be:
  #    -1/2(n/log(2*pi) + nlog(s2) + log|V| + 1/s2*((y-Xb)^t V^{-1} (y-Xb) + tau*P(b)))
  #    where tau = 2*n*lambda and |V| == 1
  #    This is maximized with respect to b when (y-Xb)^t V^{-1} (y-Xb) + tau*P(b) is minimized
  #    This is equivalent to minimizing 1/2RSS/n + lambda*P(b) which is the glmnet criteria
  #    Given \hat{b}, s2_hat = (RSS + 2*n*P(b))/n
  # The score is -2 time penalized log-Likelihood
  s2_hats = (RSSs + 2*nobs*penalties)/nobs
  scores = (nobs*log(2*pi) + nobs*log(s2_hats) + V_log_det + nobs)
  return(list(h2_index = h2_i,V_log_det = V_log_det,res.glmnet = res,scores = scores, lambda = lambda,deviance = RSSs,null_deviance = sum((y_star-mean(y_star))^2),s2 = s2_hats))
}


GridLMMnet2 = function(y,X_cov,X,V_list_setup,mc.cores,clusterType = 'mclapply',alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL,verbose = FALSE,...)
{
  nobs = length(y)
  nvars = ncol(X_cov) + ncol(X)
  if(clusterType == 'SNOW') {
    cl = makeCluster(mc.cores)
  } else{
    cl = mc.cores
  }
  
  if(is.character(V_list_setup)) {
    V_list_setup = readRDS(V_list_setup)
  } else if('setup_file' %in% names(V_list_setup)) {
    chol_V_list = V_list_setup$chol_V_files
    V_list_setup = readRDS(V_list_setup$setup)
    V_list_setup$chol_V_list = chol_V_list
  }
  h2s_matrix   = V_list_setup$h2s_matrix
  Qt           = V_list_setup$Qt
  n_SNPs_RRM   = V_list_setup$n_SNPs_RRM # used to calcualte downdate_weights based on downdate_Xs
  chol_V_list  = V_list_setup$chol_V_list
  
  penalty.factor = c(rep(0,ncol(X_cov)),rep(1,ncol(X)))
  X_full = cbind(X_cov,X)
  
  # standardize y going in - makes glmnet more stable
  mean_y = mean(y)
  sd_y = sd(y)
  y = (y-mean_y)/sd_y
  
  if(!is.null(Qt)){
    y <- as.matrix(Qt %*% y)
    X_full = Qt %*% X_full
  }
  
  registerDoParallel(cl)
  # find an appropriate range for lambda by running glmnet with nlambda=3 (minimum to get a stable value of lambda[1]) on each matrix
  if(verbose) print(sprintf('running glmnet on %d matrices',length(chol_V_list)))
  if(is.null(lambda)) {
    results_prep <- foreach::foreach(chol_V_setup = chol_V_list,.combine = rbind) %dopar% {
      run_glmnet_V(X_full,y,chol_V_setup,alpha = alpha,get_max_lambda = TRUE,penalty.factor=penalty.factor,...)
    }
    max_lambda = max(results_prep[,1])#results_prep[order(results_prep[,'score'])[1],'max_lambda']
    # max_lambda = results_prep[order(results_prep[,'score'])[1],'max_lambda']  # this seems to work
    min_lambda = max_lambda * lambda.min.ratio
    lambda = exp(seq(log(max_lambda),log(min_lambda),length=nlambda))
  }
  
  # use this lambda for all.
  results <- foreach::foreach(chol_V_setup = chol_V_list) %dopar% {
    run_glmnet_V(X_full,y,chol_V_setup,alpha = alpha,lambda=lambda,penalty.factor=penalty.factor,...)
  }
  if(inherits(cl,"cluster")) stopCluster(cl)
  
  ## select the best model by minimum score
  res = results[[1]]$res.glmnet
  h2_indexes = sapply(results,function(x) x$h2_index)
  scores = sapply(results,function(x) x$scores)
  min_score_index = apply(scores,1,function(x) which(x == min(x,na.rm=T))[1])
  min_scores = sapply(seq_along(lambda),function(i) scores[i,min_score_index[i]])
  beta = Matrix(sapply(seq_along(lambda),function(i) as.matrix(results[[min_score_index[i]]]$res.glmnet$beta[,i]))) #
  rownames(beta) = rownames(res$beta)
  deviances = sapply(seq_along(lambda),function(i) results[[min_score_index[i]]]$deviance[i])
  null_deviance = results[[1]]$null_deviance
  res$dev.ratio = 1-deviances/null_deviance
  # remove standardization of y
  res$nulldev = null_deviance * sd_y^2
  res$beta = beta * sd_y
  if(ncol(X_cov)>0 && all(X_cov[,1]==1)){
    res$a0 = res$beta[1,] + mean_y
    res$beta[1,] = 0
  }
  res$h2s = h2s_matrix[,h2_indexes[min_score_index]] # in case h2_indexes is out of order.
  res$s2s = sd_y^2 * sapply(seq_along(lambda),function(i) results[[min_score_index[i]]]$s2[i])
  res$min_scores = min_scores
  res$lambda = lambda
  res$df = colSums(res$beta != 0) - 1
  
  # clean up object. Temporary until I understand how to process these correctly
  res$npasses = NA
  
  return(res)
}