#' LASSO solutions in a linear mixed model using GridLMM
#'
#' Finds LASSO or Elastic Net solutions for a multiple regression problem with correlated errors.
#' 
#' @details Finds the full LASSO or Elastic Net solution path by running \code{\link[glmnet]{glmnet}} at each grid vertex.
#'     If \code{foldid} is provided, cross-validation scores will be calculated.
#'
#' @inheritParams GridLMM_GWAS
#' @inheritParams glmnet::glmnet
#' @param X Variables in model that well be penalized with the elastic net penalty. Covariates specified in \code{formula} are not penalized.
#' @param foldid vector of integers that divide the data into a set of non-overlapping folds for cross-validation.
#' @param ... 
#'
#' @return If \code{foldid} and \code{nfold} are null, an object with S3 class "glmnet","*" , where "*" is "elnet". See \code{\link[glmnet]{glmnet}}.
#'     Otherwise, an object with S3 class "cv.glmnet". See \code{\link[glmnet]{cv.glmnet}}.
#' @export
#'
#' @examples
GridLMMnet = function(formula,data,X, weights = NULL, 
                      centerX = TRUE,scaleX = TRUE,relmat = NULL,
                      h2_divisions = 10, h2_start = NULL,
                      alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL,
                      nfolds = NULL,foldid = NULL,
                      RE_setup = NULL, V_setup = NULL, save_V_folder = NULL,
                      diagonalize=T,svd_K = T,drop0_tol = 1e-10,mc.cores = parallel::detectCores(),clusterType = 'mclapply',verbose=T,...) 
{
  
  setup = GridLMMnet_setup(formula,data,X,weights,
                           centerX,scaleX,relmat,
                           h2_divisions,h2_start,
                           alpha,nlambda,substitute(lambda.min.ratio),lambda,
                           nfolds,foldid,
                           RE_setup,V_setup,save_V_folder,
                           diagonalize,svd_K,drop0_tol,mc.cores,clusterType,verbose,...)

  lambda.min.ratio = setup$lambda.min.ratio
  
  results = GridLMMnet2(setup,mc.cores,clusterType=clusterType,alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda=lambda,
                        prediction_weights = weights, foldid = foldid,verbose=verbose,...)
  
  # returning elnet object. Don't think it's worth including the setup object too. I can add functionality to re-use V_list if passed.
  return(results)
  
}

GridLMMnet_setup = function(formula,data,X, weights = NULL, 
                            centerX = TRUE,scaleX = TRUE,relmat = NULL,
                            h2_divisions = 10, h2_start = NULL,
                            alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL,
                            nfolds = NULL,foldid = NULL,
                            RE_setup = NULL, V_list_setup = NULL, save_V_folder = NULL,
                            diagonalize=T,svd_K = T,drop0_tol = 1e-10,mc.cores = parallel::detectCores(),clusterType = 'mclapply',verbose=T,...) {
 
  
  # -------- sort by foldid ---------- #
  if(!is.null(nfolds)) {
    if(!is.null(foldid)) {
      warning('using provided foldid, NOT nfolds')
    } else{
      foldid = sample(rep(seq(nfolds), length = nrow(data)))
    }
  }
  if(!is.null(foldid)){
    if(length(foldid) != nrow(data)) stop('wrong length of foldid')
    data$original_order = 1:nrow(data)
    new_order = order(foldid)
    data = data[new_order,,drop=FALSE]
    X = X[new_order,,drop=FALSE]
    weights = weights[new_order]
    foldid = foldid[new_order]
    diagonalize = F # can't diagonalize with cross-validation
  }
  
  
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
  intercept = 0
  if(ncol(X_cov)>0 && all(X_cov[,1]==1)) intercept = 1
  
  nvars = ncol(X_cov) + ncol(X)
  
  # -------- Random effects ---------- #
  if(is.null(RE_setup)) {
    RE_setup = make_RE_setup(formula = formula,data,relmat = relmat,verbose=verbose)
    V_setup = make_V_setup(RE_setup,weights,diagonalize,svd_K,drop0_tol,save_V_folder,verbose)
  }
  
  # -------- Form matrix of h2s to test ---------- #
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
  if(!is.null(h2_start)) {
    if(length(h2_start) != nrow(h2s_matrix)) stop('wrong length h2_start vector provided')
    if(!is.null(names(h2_start))){
      if(!all(rownames(h2s_matrix) %in% names(h2_start))) stop('missing h2s in h2_start')
      h2s_matrix = cbind(h2s_matrix,h2_start[rownames(h2s_matrix)])
    } else{
      h2s_matrix = cbind(h2s_matrix,h2_start)
    }
  }
  colnames(h2s_matrix) = NULL
  
  if(clusterType == 'mclapply') {
    registerDoParallel(mc.cores)
  } else{
    registerDoParallel(1)
    if(verbose) pb = txtProgressBar(min=0,max = length(chol_V_list),style=3)
  }
  V_setup$chol_V_list = foreach(h2s = iter(h2s_matrix,by='col')) %dopar% {
    h2s = h2s[1,]
    chol_V_setup = make_chol_V_setup(V_setup,h2s)
    if(!is.null(V_setup$save_V_folder)) chol_V_setup = chol_V_setup$file  # only store file name if save_V_folder is provided
    if(verbose && exists('pb')) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    return(chol_V_setup)
  }
  if(verbose && mc.cores == 1 && exists('pb')) close(pb)
  
  penalty.factor = c(rep(0,ncol(X_cov)),rep(1,ncol(X)))
  X_full = cbind(X_cov,X)
  
  # standardize y going in - makes glmnet more stable
  mean_y = mean(y)
  sd_y = sd(y)
  y = (y-mean_y)/sd_y
  
  Qt = V_setup$Qt
  if(!is.null(Qt)){
    y <- as.matrix(Qt %*% y)
    X_full = as.matrix(Qt %*% X_full)
  }
  
  setup = list(
    y = y,
    X_full = X_full,
    intercept = intercept,
    penalty.factor = penalty.factor,
    foldid = foldid,
    data = data,
    mean_y = mean_y,
    sd_y = sd_y,
    weights = weights,
    lambda.min.ratio = eval(lambda.min.ratio),
    V_setup = V_setup,
    h2s_matrix = h2s_matrix
  )
  
  return(setup)  
   
}


run_glmnet_V = function(X_full,y, chol_V, alpha = alpha, 
                        lambda = NULL,nlambda, lambda.min.ratio,get_max_lambda = FALSE,penalty.factor,holdOut = NULL,...){
  nobs = length(y)
  
  if(!is.null(holdOut)){
    if(length(holdOut) != nobs) stop('wrong length of foldid')
    if(sum(diff(holdOut)<0) > 1) stop('foldid must be sorted')
    
    # first pull out validation set
    y_sub = y[holdOut]
    X_sub = X_full[holdOut,,drop=FALSE]
    chol_V_full = chol_V
    
    # now reduce data to only training set to fit model
    y = y[!holdOut]
    X_full = X_full[!holdOut,,drop=FALSE]
    chol_V = t(chol_dropRows(t(as.matrix(chol_V)),min(which(holdOut)),sum(holdOut)))
  }
  
  V_log_det = 2*sum(log(diag(chol_V)))
  
  chol_V_inv = as.matrix(solve(chol_V))
  
  # calculate a normalization constant based on V_inv
  ci = sum(chol_V_inv^2)/nobs # normalization constant based on trace of V_inv
  
  # normalize chol_V_inv
  chol_V_inv = chol_V_inv/sqrt(ci)
  V_log_det = -2*sum(log(diag(chol_V_inv)))
  
  y_star = crossprod_cholR(chol_V_inv,y)
  X_star = crossprod_cholR(chol_V_inv,X_full)
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
  errors = y_star[,1] - as.matrix(X_star %*% res$beta)
  RSSs = colSums(errors^2)
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
  
  # predicting in held-out set
  prediction_errors = NULL
  if(!is.null(holdOut)){
    # if(sum(chol_V != 0) > nobs) recover()
    # want X_1*beta_hat + V_12 * V_22^{-1}(y-X_2*b_hat)
    # note: chol_V_full %*% t(chol_V_full) * ci = V_22^{-1}
    
    predicted = as.matrix(X_sub %*% res$beta) + as.matrix((t(chol_V_full[,holdOut]) %*% chol_V_full[,!holdOut]) %*% chol_V_inv %*% errors) * ci
    prediction_errors = y_sub - predicted
  }
  return(list(V_log_det = V_log_det,res.glmnet = res,scores = scores, lambda = lambda,
              deviance = RSSs,null_deviance = sum((y_star-mean(y_star))^2),s2 = s2_hats,
              prediction_errors = prediction_errors))
}


get_lambda_sequence = function(setup,alpha = 1,nlambda = 100,lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001),verbose = TRUE,...) {
  
  y = setup$y
  X_full = setup$X_full
  V_setup = setup$V_setup
  penalty.factor = setup$penalty.factor
  chol_V_list  = V_setup$chol_V_list
  nobs = length(y)
  nvars = ncol(X_full)
  
  
  # find an appropriate range for lambda by running glmnet with nlambda=3 (minimum to get a stable value of lambda[1]) on each matrix
  if(verbose) print(sprintf('calculating lambda sequence for %d matrices',length(chol_V_list)))
  results_prep <- foreach::foreach(chol_V_setup = chol_V_list,.combine = rbind,.packages = c('glmnet')) %dopar% {
    if(is.character(chol_V_setup)) {
      chol_V_setup <- readRDS(chol_V_setup)
    }
    run_glmnet_V(X_full,y,chol_V_setup$chol_V,alpha = alpha,get_max_lambda = TRUE,penalty.factor=penalty.factor,...)
  }
  max_lambda = max(results_prep[,1])#results_prep[order(results_prep[,'score'])[1],'max_lambda']
  # max_lambda = results_prep[order(results_prep[,'score'])[1],'max_lambda']  # this seems to work
  min_lambda = max_lambda * lambda.min.ratio
  lambda = exp(seq(log(max_lambda),log(min_lambda),length=nlambda))
  
  return(lambda)
}


run_GridLMMnet = function(setup,alpha = 1,lambda,verbose = TRUE,...) {
  y = setup$y
  X_full = setup$X_full
  foldid = setup$foldid
  V_setup = setup$V_setup
  penalty.factor = setup$penalty.factor
  chol_V_list  = V_setup$chol_V_list
  
  if(verbose) {
    print(sprintf('running glmnet on %d grid cells',length(chol_V_list)))
    pb = txtProgressBar(min=0,max = length(chol_V_list),style=3)
  }
  results <- foreach::foreach(chol_V_setup = chol_V_list, h2_index = seq_along(chol_V_list),.packages = c('glmnet')) %dopar% {
    if(is.character(chol_V_setup)) {
      chol_V_setup <- readRDS(chol_V_setup)
    }
    res = run_glmnet_V(X_full,y,chol_V_setup$chol_V,alpha = alpha,lambda=lambda,penalty.factor=penalty.factor,...)
    res$h2_index = h2_index
    if(!is.null(foldid)){
      nfolds = max(foldid)
      res_list = foreach::foreach(i = seq(nfolds)) %do% {
        holdOut = foldid == i
        res_i = run_glmnet_V(X_full,y,chol_V_setup$chol_V,alpha = alpha,lambda=lambda,penalty.factor=penalty.factor,holdOut = holdOut,...)
        return(res_i[c('scores','prediction_errors')])
      }
      res$cv_scores = sapply(res_list,function(x) x$scores)  # score for each fold for each lambda
      res$prediction_errors = lapply(res_list,function(x) x$prediction_errors)
    }
    
    if(verbose) {
      setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    }
    
    res
  }
  if(verbose) close(pb)
  return(results)
}

collect_results_GridLMMnet = function(setup,results,lambda) {
  mean_y = setup$mean_y
  sd_y = setup$sd_y
  h2s_matrix = setup$h2s_matrix
  foldid = setup$foldid
  prediction_weights = setup$weights
  nobs = length(setup$y)
  
  res = results[[1]]$res.glmnet
  h2_indexes = sapply(results,function(x) x$h2_index)
  
  ## select the best model by minimum score
  scores = sapply(results,function(x) x$scores)
  min_score_index = apply(scores,1,which.min)
  min_scores = sapply(seq_along(lambda),function(i) scores[i,min_score_index[i]])
  beta = Matrix(sapply(seq_along(lambda),function(i) as.matrix(results[[min_score_index[i]]]$res.glmnet$beta[,i]))) #
  rownames(beta) = rownames(res$beta)
  deviances = sapply(seq_along(lambda),function(i) results[[min_score_index[i]]]$deviance[i])
  null_deviance = results[[1]]$null_deviance
  res$dev.ratio = 1-deviances/null_deviance
  # remove standardization of y
  res$nulldev = null_deviance * sd_y^2
  res$beta = beta * sd_y
  if(setup$intercept == 1){
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
  
  
  if(!is.null(foldid)){
    # this is duplicating cv.glmnet, so returns a differently structured object
    # for each lambda, then for each fold, select best h2 by score, pull out this cvm
    # select a lambda
    min_scores_index_cv = sapply(seq_along(lambda),function(l) {
      # go through all h2's, pull out score vectors for each fold: (n(fold) x n(h2) matrix)
      lambda_scores_cv = sapply(results,function(x) x$cv_scores[l,])
      # for each fold, find the index of the smallest (best) score
      apply(lambda_scores_cv,1,which.min)
    })
    # pull out prediction_errors for selected models for each fold
    cvraw = sapply(seq_along(lambda),function(l) {
      indexes = min_scores_index_cv[,l]
      do.call(c,lapply(seq_along(indexes),function(i) results[[indexes[[i]]]]$prediction_errors[[i]][,l]))^2
    })
    # remove standardization of y
    cvraw = cvraw*sd_y^2
    # weight cvraw
    if(is.null(prediction_weights)) prediction_weights = rep(1,nobs)
    cvm = apply(cvraw, 2, weighted.mean, w = prediction_weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                      w = prediction_weights, na.rm = TRUE)/(nobs - 1))
    res = list(
      lambda = lambda,
      cvm = cvm,
      cvsd = cvsd,
      cvup = cvm + cvsd,
      cvlo = cvm - cvsd,
      nzero = colSums(res$beta != 0),
      name = "Mean-Squared Error",
      glmnet.fit = res,
      foldid = foldid
    )
    res = c(res,glmnet::getmin(lambda,res$cvm,res$cvsd))
    class(res) = 'cv.glmnet'
  }
  return(res)
}

start_cluster = function(mc.cores = my_detectCores(),type = 'mclapply',...) {
  if(type == 'mclapply') {
    cl = mc.cores
  } else{
    cl = makeCluster(mc.cores,type = type,...)
  }
  registerDoParallel(cl)
  cl
}

GridLMMnet2 = function(setup,mc.cores,clusterType = 'mclapply',alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL,
                       prediction_weights = NULL,foldid = NULL,verbose = FALSE,...)
{
  
  y = setup$y
  X_full = setup$X_full
  V_setup = setup$V_setup
  foldid = setup$foldid
  penalty.factor = setup$penalty.factor
  mean_y = setup$mean_y
  sd_y = setup$sd_y
  prediction_weights = setup$weights
  
  h2s_matrix   = V_setup$h2s_matrix
  chol_V_list  = V_setup$chol_V_list
  
  nobs = length(y)
  
  # find lambda sequence
  if(is.null(lambda)) {
    cl = start_cluster(mc.cores,clusterType)
    lambda = get_lambda_sequence(setup,alpha = alpha,nlambda = nlambda,lambda.min.ratio = lambda.min.ratio,verbose = verbose,...)
    if(inherits(cl,"cluster")) stopCluster(cl)
  }
  
  # run glmnet across the grid, including cross-validation
  cl = start_cluster(mc.cores,clusterType)
  results = run_GridLMMnet(setup,alpha = alpha,lambda = lambda,verbose = verbose,...)
  if(inherits(cl,"cluster")) stopCluster(cl)
  
  # collect all results
  res = collect_results_GridLMMnet(setup,results,lambda)
  
  return(res)
}
