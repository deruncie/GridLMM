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
GridLMMnet = function(formula,data,X, X_ID = 'ID', weights = NULL, 
                      centerX = TRUE,scaleX = TRUE,relmat = NULL,
                      h2_step = 0.1, h2_start = NULL,
                      alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL,
                      penalty.factor = NULL,
                      nfolds = NULL,foldid = NULL,
                      RE_setup = NULL, V_setup = NULL, save_V_folder = NULL,
                      diagonalize=T,mc.cores = parallel::detectCores(),clusterType = 'mclapply',verbose=T,...) 
{
  
  # ----------- setup model ------------- #
  setup = GridLMMnet_setup(formula,data,X,X_ID, weights,
                           centerX,scaleX,relmat,
                           alpha,nlambda,substitute(lambda.min.ratio),lambda,
                           penalty.factor,
                           nfolds,foldid,
                           RE_setup,V_setup,save_V_folder,
                           diagonalize,mc.cores,clusterType,verbose,...)
  nobs = length(setup$y)
  nvars = ncol(setup$X_full)
  # recover()
  
  # ----------- setup Grid ------------- #
  setup$h2s_matrix = setup_Grid(names(setup$V_setup$RE_setup),h2_step,h2_start)
  
  cl = start_cluster(mc.cores,clusterType)
  setup$V_setup = calculate_Grid(setup$V_setup,setup$h2s_matrix,verbose,invIpKinv = T)
  if(inherits(cl,"cluster")) stopCluster(cl)

  # ----------- run model ------------- #

  # find lambda sequence
  lambda = setup$lambda
  if(is.null(lambda)) {
    cl = start_cluster(mc.cores,clusterType)
    lambda = get_lambda_sequence(setup,alpha = alpha,nlambda = nlambda,lambda.min.ratio = lambda.min.ratio,verbose = verbose,...)
    if(inherits(cl,"cluster")) stopCluster(cl)
  }
  
  # run glmnet across the grid, including cross-validation
  cl = start_cluster(mc.cores,clusterType)
  grid_results_list = run_GridLMMnet(setup,alpha = alpha,lambda = lambda,verbose = verbose,...)
  if(inherits(cl,"cluster")) stopCluster(cl)
  
  # collect all results
  res = collect_results_GridLMMnet(setup,grid_results_list,lambda)
  res$lambda = res$lambda * setup$sd_y^2
  
  return(res)
  
  
}

GridLMMnet_setup = function(formula,data,X, X_ID = 'ID', weights = NULL, 
                            centerX = TRUE,scaleX = TRUE,relmat = NULL,
                            alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL,
                            penalty.factor = NULL,
                            nfolds = NULL,foldid = NULL,
                            RE_setup = NULL, V_setup = NULL, save_V_folder = NULL,
                            diagonalize=T,mc.cores = parallel::detectCores(),clusterType = 'mclapply',verbose=T,...) {
 
  
  # -------- sort by foldid ---------- #
  original_order = NULL
  if(!is.null(nfolds)) {
    if(!is.null(foldid)) {
      warning('using provided foldid, NOT nfolds')
    } else{
      foldid = sample(rep(seq(nfolds), length = nrow(data)))
    }
  }
  if(!is.null(foldid)){
    if(length(foldid) != nrow(data)) stop('wrong length of foldid')
    original_order = 1:nrow(data)  # NOT SAFE!!!
    new_order = order(foldid)
    
    # check that all variables are in data. If not, add
    for(term in all.vars(formula)) {
      if(!term %in% colnames(data)) {
        data[[term]] = eval(parse(text=term))
      }
    }
    data = data[new_order,,drop=FALSE]
    X = X[new_order,,drop=FALSE]
    if(!is.null(weights)) weights = weights[new_order]
    foldid = foldid[new_order]
    original_order = original_order[new_order]
    diagonalize = F # can't diagonalize with cross-validation
  }
  # -------- prep Mixed Models ---------- #
  MM = prepMM(formula,data,weights,other_formulas = NULL,
              relmat,X,X_ID,proximal_markers=NULL,V_setup,diagonalize, svd_K = TRUE,drop0_tol = 1e-10,save_V_folder, verbose)
  lmod = MM$lmod
  RE_setup = MM$RE_setup
  V_setup = MM$V_setup
  
  y = matrix(lmod$fr[,1])
  n = nobs = nrow(y)
  X_cov = lmod$X
  intercept = 0
  if(ncol(X_cov)>0 && all(X_cov[,1]==1)) {
    intercept = 1
  } else if(ncol(X_cov) == 0 && all(X[,1] == 1)) {
    intercept = 1
  }
  data = lmod$fr
  
  Z_X = model.matrix(formula(sprintf("~0+%s",X_ID)),droplevels(data))
  colnames(Z_X) = sub(X_ID,'',colnames(Z_X),fixed = T)
  if(is.null(rownames(X))) stop(sprintf('X must have rownames that correspond with column %s in data',X_ID))
  stopifnot(all(colnames(Z_X) %in% rownames(X)))
  X = Z_X %*% X[colnames(Z_X),]
  
  nvars = ncol(X_cov) + ncol(X)
  if(!is.null(penalty.factor)) {
    if(length(penalty.factor) != ncol(X)) stop("Wrong length of penalty.factor")
    penalty.factor = c(rep(0,ncol(X_cov)),penalty.factor)
  } else{
    penalty.factor = c(rep(0,ncol(X_cov)),rep(1,ncol(X)))
  }
  X_full = cbind(X_cov,X)
  
  # standardize y going in - makes glmnet more stable
  mean_y = mean(y)
  sd_y = 1#sd(y)*sqrt((n-1)/n)
  y = (y-mean_y)/sd_y
  # if(!is.null(lambda)) lambda = lambda * sd_y^2
  
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
    lambda = lambda,
    V_setup = V_setup,
    original_order = original_order
  )
  
  return(setup)  
   
}


make_chol_V_setup = function(V_setup,h2s,invIpKinv = F){
  if(!is.null(V_setup$save_V_folder)) {
    chol_V_files = read.csv(sprintf('%s/chol_V_index.csv',V_setup$save_V_folder),stringsAsFactors = FALSE)
    if(nrow(chol_V_files)>0){
      # identify the chol_V_setup file with match h2s and load
      h2_matches = rowSums(abs(sweep(chol_V_files[,-1,drop=FALSE],2,h2s,'-')))
      h2_match = which(h2_matches<1e-10)
      if(length(h2_match) == 1){
        chol_V_setup = readRDS(chol_V_files$file[h2_match])
        return(chol_V_setup)
      }
    }
  }
  downdate_ratios = V_setup$downdate_ratios
  total_K = 0
  for(i in 1:length(h2s)) total_K = total_K + h2s[i] * as.matrix(V_setup$ZKZts[[i]]) * downdate_ratios[i]
  V = total_K + (1-sum(h2s)) * V_setup$resid_V
  # test if V is diagonal or sparse
  non_zero_upper_tri_V = abs(V[upper.tri(V,diag = FALSE)]) > 1e-10
  if(sum(non_zero_upper_tri_V) == 0) {  # V is diagonal
    chol_V = as(diag(sqrt(diag(V))),'CsparseMatrix')
  } else if(sum(non_zero_upper_tri_V) < length(non_zero_upper_tri_V) / 2) { # V is sparse
    V = Matrix(V)
    chol_V = as(chol(V),'CsparseMatrix')
  } else{ # V is dense, use Eigen LLT
    chol_V = chol_c(V)
  }
  if(inherits(chol_V,'dtCMatrix')) chol_V = as(chol_V,'dgCMatrix')
  V_log_det = 2*sum(log(diag(chol_V)))
  
  chol_V_setup = list(chol_V = chol_V,V_log_det = V_log_det, h2s = h2s)
  
  if(invIpKinv) {
    # want (I + e/aK^{-1})^{-1}
    # a/e*K - a/e*K(I+a/e*K)^{-1}a/e*K
    # V = e/(a+e)I+a/(a+e)K
    # (a+e)/eV = I+a/e*K
    # a/e*K - a/e*K(e/(a+e)*V^{-1})a/e*K
    # a/e*(K - e/(a+e)*(a/e)*crossprod(solve(t(R),K)))
    a = sum(h2s)
    e = 1-a
    h2 = sum(h2s)
    if(h2 == 0) {
      total_K = 0*total_K
    } else{ 
      total_K = total_K/h2
    }
    chol_V_setup$invIpKinv = (h2/(1-h2)*(total_K - h2* crossprod(solve(t(chol_V),total_K))))
    # sK = svd(K)
    if(sum(is.infinite(total_K))>0) recover()
    # chol_V_setup$sK <- svd(total_K)
    # chol_K = suppressWarnings(chol(total_K,pivot=T))
    if(h2 == 0) {
      chol_V_setup$L = matrix(0,nrow = nrow(total_K),ncol = 1)
    } else{
      ldlK = BSFG::LDLt(total_K)
      d = sum(ldlK$d > 1e-10)
      chol_V_setup$L = ldlK$P %*% ldlK$L[,1:d] %*% diag(sqrt(ldlK$d[1:d]))
    }
  }
  
  if(!is.null(V_setup$save_V_folder)) {
    # save chol_V_setup as a file, append the file to the chol_V_list csv
    chol_V_file = sprintf('%s/chol_V_%s.rds',V_setup$save_V_folder,paste(h2s,collapse='_'))
    saveRDS(chol_V_setup,file = chol_V_file)
    system(sprintf('echo "%s,%s" >> %s/chol_V_index.csv',chol_V_file,paste(h2s,collapse=','),V_setup$save_V_folder))
    chol_V_setup$file = chol_V_file
  }
  return(chol_V_setup)
}


calculate_Grid = function(V_setup,h2s_matrix,verbose = TRUE,invIpKinv = F){
  
  if(verbose) {
    sprintf('Generating V decompositions for %d grid cells', ncol(h2s_matrix))
    pb = txtProgressBar(min=0,max = ncol(h2s_matrix),style=3)
  }
  
  V_setup$chol_V_list = foreach(h2s = iter(t(h2s_matrix),by='row')) %dopar% {
    h2s = h2s[1,]
    chol_V_setup = make_chol_V_setup(V_setup,h2s,invIpKinv)
    if(!is.null(V_setup$save_V_folder)) chol_V_setup = chol_V_setup$file  # only store file name if save_V_folder is provided
    if(verbose && exists('pb')) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    return(chol_V_setup)
  }
  if(verbose && exists('pb')) close(pb)
  
  return(V_setup)
}


run_glmnet_V = function(X_full,y, chol_V_setup, alpha = alpha, 
                        lambda = NULL,nlambda, lambda.min.ratio,get_max_lambda = FALSE,penalty.factor,holdOut = NULL,...){
  nobs = length(y)
  h2 = sum(chol_V_setup$h2s)
  chol_V = chol_V_setup$chol_V / sqrt(1-h2)
  
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
  
  # penalty = tr_Vinv / nobs * lambda
  tr_Vinv = sum(chol_V_inv^2)
  
  y_star = crossprod_cholR(chol_V_inv,y)
  X_star = crossprod_cholR(chol_V_inv,X_full)
  sd_y_star = sd(y_star)*sqrt((nobs-1)/nobs)
  
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
    resid_y_star = y_star - mean(y_star)
    if(any(penalty.factor == 0)) {
      resid_y_star = resid(lm(y_star~X_star[,penalty.factor == 0]))
    }
    s2_hat = sum(resid_y_star^2)/nobs
    
    # Note: if alpha == 0, this fails. Set alpha = max(alpha,0.01)
    alpha = max(alpha,0.001)

    # max_lambda = nobs/tr_Vinv*max(penalty.factor*sum(penalty.factor)/(length(penalty.factor)-1)*abs(t(resid_y_star) %*% X_star ) / ( alpha * nobs))
    max_lambda = max(penalty.factor*sum(penalty.factor)/(length(penalty.factor)-1)*abs(t(resid_y_star) %*% X_star ) / ( alpha * nobs))
    # max_lambda = 1/(1-h2)*max(penalty.factor*sum(penalty.factor)/(length(penalty.factor)-1)*abs(t(resid_y_star) %*% X_star ) / ( alpha * nobs))
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
    # lambda = nobs / tr_Vinv * res$lambda * sd_y_star^2
    lambda = res$lambda * sd_y_star^2
    # alternatively, use family = 'mgaussian' and standarize.response = F
    # res = glmnet(X_star,y_star,family = 'mgaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,standardize.response = F,
    #              alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,intercept = intercept,standardize = FALSE,...)
    # class(res)[1] = 'elnet' # so we can use plot.elnet
    # lambda = res$lambda
  } else{
    # lambda = lambda * ci
    # res = glmnet(X_star/sd_y_star,y_star/sd_y_star,family = 'gaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,
    #              alpha = alpha, lambda = tr_Vinv / nobs * lambda/sd_y_star^2 ,intercept = intercept,standardize = FALSE,...)
    res = glmnet(X_star/sd_y_star,y_star/sd_y_star,family = 'gaussian',weights = rep(1,length(y_star)),#penalty.factor = penalty.factor,
                 alpha = alpha, lambda = lambda/sd_y_star^2,intercept = intercept,standardize = FALSE,...)
    # res = glmnet(X_star/sd_y_star,y_star/sd_y_star,family = 'gaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,
    #              alpha = alpha, lambda = lambda/sd_y_star^2,intercept = intercept,standardize = FALSE,...)
    # res = glmnet(X_star/sd_y_star,y_star/sd_y_star,family = 'gaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,
    #              alpha = alpha, lambda = lambda/sd_y_star^2*(1-h2) ,intercept = intercept,standardize = FALSE,...)
    # alternatively, use family = 'mgaussian' and standarize.response = F
    # res = glmnet(X_star,y_star,family = 'mgaussian',weights = rep(1,length(y_star)),penalty.factor = penalty.factor,standardize.response = F,
    #               alpha = alpha, lambda = lambda,intercept = intercept,standardize = FALSE,...)
    # class(res)[1] = 'elnet' # so we can use plot.elnet
  }
  if(intercept){
    res$beta[1,] = res$a0*sd_y_star
    res$a0[] = 0
  }
  # res$lambda = lambda
  # plot(res,'lambda',main = h2)
  errors = y_star[,1] - as.matrix(X_star %*% res$beta)
  RSSs = colSums(errors^2)
  # recover()
  penalty.factor = penalty.factor*(length(penalty.factor)-1)/sum(penalty.factor) # adjust penalty.factor to sum to p-1 (for intercept)
  # penalty.factor = penalty.factor*(length(penalty.factor))/sum(penalty.factor) # adjust penalty.factor to sum to p-1 (for intercept)
  penalties = colSums(penalty.factor*((1-alpha)/2*res$beta^2 + alpha*abs(res$beta)))# from: https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
  
  # V = a/(a+e)*K + e/(a+e)I
  # (a+e)V = aK+eI
  # V2 = a/e K + I
  # eV2 = aK + eI
  # V2 = (a+e)/e * V = 1/(1-h2)*V
  
  dets=0
  # recover()
  # s2_hats = (RSSs + 2*tr_Vinv*lambda*penalties)/nobs
  # s2_hats = (RSSs + 2*nobs*lambda*penalties)/nobs
  # scores = (nobs*log(2*pi) + nobs*log(s2_hats) + V_log_det + nobs)
  
  s2_hats = (RSSs + 2*nobs*lambda*penalties)/(nobs+2*p)
  scores = (nobs+2*p)*log(s2_hats) + V_log_det
  
  #GCV
  # 1/nRSS/(1-q/n)^2
  #q = tr(X(XtX+lambdaW-)^{-1}Xt)
  # recover()
  # qt = sapply(seq_along(lambda),function(i) {
  #   XWXt = 1/lambda[i]*X_star %*% (abs(res$beta[,i])*t(X_star))
  #   H = XWXt - XWXt %*% solve(XWXt+diag(1,nobs),XWXt)
  #   sum(diag(H))
  # })
  # k = colSums(res$beta != 0)
  # scores = 1/nobs*RSSs/(1-qt/n)^2
  
  # scX = svd(X_star)
  # k = pmin(colSums(res$beta != 0),nobs)
  # k = ncol(X_star)
  # s2_hats = (RSSs + 2*nobs*lambda*penalties)/(nobs-k)
  # dets = 2*sum(log(scX$d[scX$d > 1e-10]))
  # scores = (nobs-k)*log(s2_hats) + V_log_det + nobs - k + dets
  
  # s2_hats = ((1-h2)*RSSs + 2*nobs*lambda*penalties)/nobs
  # scores = (nobs*log(2*pi) + nobs*log(s2_hats) + V_log_det + nobs - nobs*log(1-h2))
  # cXs = crossprod(X_star)
  # tcXs = tcrossprod(X_full)
  # V = crossprod(chol_V)
  # dets = sapply(lambda,function(l)  {
  #   if(nobs > p) {
  #     -V_log_det + determinant(V + tcXs/(nobs*l))$modulus
  #   } else {
  #     determinant(diag(1,ncol(X_full)) + cXs/(nobs*l))$modulus
  #   }
  #   })
  # dets=0
  # dets=0
  
  # # The following seems to be right for ridge, but not lasso.
  # dets = 0
  # if(alpha < 1) {
  #   scX = svd(X_star)
  #   d = scX$d
  #   if(length(d) < p) d = c(d,rep(0,(p-length(d))))
  #   dets = dets + (1-alpha)*sapply(lambda,function(l) {
  #     # determinant(diag(1,p) + diag(d^2/(nobs*l)))$modulus
  #     # sum(log(1+d^2/(tr_Vinv*l)))
  #     sum(log(1+d^2/(nobs*l)))
  #   })
  # }
  # if(alpha > 0) {
  #   dets = dets + alpha * sapply(seq_along(lambda),function(i) {
  #     coefs = which(res$beta[,i] != 0)
  #     s_Xi = svd(X_full[,coefs])
  #     s_Xstari = svd(X_star[,coefs])
  #     2*sum(log(s_Xstari$d))# - 2*sum(log(s_Xi$d))
  #   })
  # }
  # # s2_hats = (RSSs + 2*tr_Vinv*lambda*penalties)/nobs
  # s2_hats = (RSSs + 2*nobs*lambda*penalties)/nobs
  # scores = (nobs*log(2*pi) + nobs*log(s2_hats) + V_log_det + nobs + dets)
  # 
  
  # The following is trying to derive it from the posterior
  # XtX = crossprod(X_star)
  # dets1 = sapply(lambda,function(l) {
  #     XtXi = XtX
  #     diag(XtXi) = diag(XtXi) + 2*nobs*l
  #     determinant(XtXi)$mod
  # })
  # d = scX$d
  # if(length(d) < p) d = c(d,rep(0,(p-length(d))))
  # # if(length(d) < n) d = c(d,rep(0,(n-length(d))))
  # dets = sapply(lambda,function(l) {
  #   sum(log(d^2 + nobs*l))
  # })
  # Vtb = t(scX$v) %*% res$beta
  # # mod_penalties1 = sapply(seq_along(lambda),function(i) {
  # #   XtXi = XtX
  # #   diag(XtXi) = diag(XtXi) + 2*nobs*lambda[i]
  # #   matrix(res$beta[,i],nr=1) %*% XtXi %*% res$beta[,i]
  # # })
  # mod_penalties = sapply(seq_along(lambda),function(i) {
  #   # XtXi = XtX
  #   # diag(XtXi) = diag(XtXi) + 2*nobs*lambda[i]
  #   # matrix(res$beta[,i],nr=1) %*% solve(XtXi,res$beta[,i])
  #   sum(((sqrt(scX$d^2 + nobs*lambda[i]))*Vtb[,i])^2)
  # })
  # # recover()
  # s2_hats = (sum(y_star^2) - mod_penalties)/nobs
  # scores = nobs*log(2*pi) + nobs*log(s2_hats) + nobs + V_log_det + dets - p*log(n*lambda)
  # s2_hats = (RSSs + 2*nobs*lambda*penalties)/nobs
  # s2_hats = (RSSs + 2*nobs*lambda*penalties + mod_penalties)/(p)
  # scores = ((p)*log(2*pi) + (p)*log(s2_hats) + nobs +dets) #+ V_log_det

  # try using maximum of random effects
  # recover()
  # e = y[,1] - as.matrix(X_full %*% res$beta)
  # s2_e = s2_hats#*(1-sum(chol_V_setup$h2s))
  # s2_a = s2_hats*h2/(1-h2)
  # s2_e = s2_hats*(1-h2)
  # s2_a = s2_hats*h2
  
   # recover()
  # scores = n*log(s2_e) + n*log(s2_a) +V_log_det + (colSums(e^2) - diag(t(e) %*% chol_V_setup$invIpKinv %*% e))/(s2_e) + 2*nobs*lambda*penalties/s2_hats
  # scores = nobs*log(s2_e) + nobs*log(s2_a) + colSums(e^2)/s2_e - diag(t(e) %*% chol_V_setup$invIpKinv %*% e)/s2_e + 2*nobs*lambda*penalties/s2_hats
  # recover()
  # try({
  # b2 = solve(t(X_full) %*% X_full + diag(s2_e/s2_a,ncol(X_full)),t(X_full) %*% e)
  # res$beta = res$beta + b2
  # })
  # d = chol_V_setup$sK$d
  # k = sum(d > 1e-10)
  # d = d[1:k]
  # U = chol_V_setup$sK$u[,1:k]
  # Z = t(sqrt(d)*t(U))
  # b = 1/(diag(crossprod(Z)) + (1-h2)/h2) * t(Z) %*% e 
  # e2 = e - Z %*% b
  # scores = nobs*log(s2_e) + k*log(s2_a) + colSums(e2)/s2_e + colSums(b^2)/s2_a
  
  # if(h2 > 0) {
  #   k = ncol(chol_V_setup$L)
  #   b = solve(crossprod(chol_V_setup$L) + diag((1-h2)/h2,ncol(chol_V_setup$L)),crossprod(chol_V_setup$L,e))
  #   e2 = e - chol_V_setup$L %*% b
  #   scores = nobs*log(s2_e)+ k*log(s2_a) + (colSums(e2^2) + 2*nobs*lambda*penalties)/s2_e + colSums(b^2)/s2_a
  # } else {
  #   scores = nobs*log(s2_e) + (colSums(e^2) + 2*nobs*lambda*penalties)/s2_e
  # }
  
  # yty-2yta+ata + atKinva
  # yty - 2ytWy + at(I+Kinv)a
  # yty - 2ytWy + ytWy
  # 1/s2b = 1/s2*(s2/s2b) = lambda/s2
  # 1/eI + 1/aKinv
  # 1/e(I+e/aKinv)
  # a_hat = (I+e/a*Kinv)^{-1}e
  #the following is an attempt for lasso, but doesn't seem to work
  # dets = sapply(seq_along(lambda),function(i)  {
  #   j = res$beta[,i] != 0 & penalty.factor > 0
  #   determinant(diag(1,sum(j)) + crossprod(X_star[,j])/(tr_Vinv*lambda[i]))$modulus
  # })
  # 
  # scores = (nobs*log(2*pi) + nobs*log(s2_hats) + V_log_det + nobs + dets)
  # # scores = ((nobs+p)*log(2*pi) + (nobs+p)*log(s2_hats) + V_log_det + nobs - p*log(nobs/(lambda*tr_Vinv)))
  # s2_hats = (RSSs + 2*nobs*lambda*penalties)/(nobs+p)
  # scores = ((nobs+p)*log(2*pi) + (nobs+p)*log(s2_hats) + V_log_det + nobs - p*log(1/(2*nobs*lambda)))
# recover()
  # predicting in held-out set
  prediction_errors = NULL
  if(!is.null(holdOut)){
    # if(sum(chol_V != 0) > nobs) recover()
    # want X_1*beta_hat + V_12 * V_22^{-1}(y-X_2*b_hat)
    # note: chol_V_full %*% t(chol_V_full) * ci = V_22^{-1}
    
    predicted = as.matrix(X_sub %*% res$beta) + as.matrix((t(chol_V_full[,holdOut]) %*% chol_V_full[,!holdOut]) %*% chol_V_inv %*% errors)# * ci
    prediction_errors = y_sub - predicted
  }
  return(list(V_log_det = V_log_det,res.glmnet = res,scores = scores, lambda = lambda,
              deviance = RSSs,null_deviance = sum((y_star-mean(y_star))^2),s2 = s2_hats,dets = dets,
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
  if(verbose) print(sprintf('calculating lambda sequence for %d grid cells',length(chol_V_list)))
  results_prep <- foreach::foreach(chol_V_setup = chol_V_list,.combine = rbind,.packages = c('glmnet')) %dopar% {
    if(is.character(chol_V_setup)) {
      chol_V_setup <- readRDS(chol_V_setup)
    }
    run_glmnet_V(X_full,y,chol_V_setup,alpha = alpha,get_max_lambda = TRUE,penalty.factor=penalty.factor,...)
  }
  results_prep = as.matrix(results_prep,nrow = length(chol_V_list))
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
    res = run_glmnet_V(X_full,y,chol_V_setup,alpha = alpha,lambda=lambda,penalty.factor=penalty.factor,...)
    res$h2_index = h2_index
    if(!is.null(foldid)){
      nfolds = max(foldid)
      res_list = foreach::foreach(i = seq(nfolds)) %do% {
        holdOut = foldid == i
        res_i = run_glmnet_V(X_full,y,chol_V_setup,alpha = alpha,lambda=lambda,penalty.factor=penalty.factor,holdOut = holdOut,...)
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
  # recover()
  if(is.null(foldid)) {
    scores = sapply(results,function(x) x$scores)
    min_score_index = apply(scores,1,which.min)
    min_scores = sapply(seq_along(lambda),function(i) scores[i,min_score_index[i]])
  } else{
    # recover()
    # pull out prediction_errors for selected models for each fold
    
    cv_results = lapply(1,function(l) {
      cvraw_lambda = sapply(seq_along(h2_indexes),function(i) {
        do.call(c,lapply(seq_along(unique(foldid)),function(f) results[[i]]$prediction_errors[[f]][,l]^2))
      })
      # remove standardization of y
      cvraw_lambda = cvraw_lambda*sd_y^2
      # weight cvraw
      if(is.null(prediction_weights)) prediction_weights = rep(1,nobs)
      # group by foldid
      cvob = glmnet::cvcompute(cvraw_lambda,prediction_weights,foldid,rep(length(h2_indexes),max(foldid)))
      cvraw_lambda = cvob$cvraw
      weights = cvob$weights
      N = cvob$N
      cvm_lambda = apply(cvraw_lambda, 2, weighted.mean, w = weights, na.rm = TRUE)
      cvsd_lambda = sqrt(apply(scale(cvraw_lambda, cvm_lambda, FALSE)^2, 2, weighted.mean,
                        w = weights, na.rm = TRUE)/(N - 1))
      list(cvm_lambda=cvm_lambda,cvsd_lambda=cvsd_lambda,cvraw_lambda=cvraw_lambda)
    })
    cvms = sapply(cv_results,function(x) x$cvm_lambda)
    min_score_index = sapply(cv_results,function(x) which.min(x$cvm_lambda))
    cvraw_lambda_base = cv_results[[1]]$cvraw_lambda[,1]
    cv_results = lapply(seq_along(lambda),function(l) {
      cvraw_lambda = sapply(seq_along(h2_indexes),function(i) {
        do.call(c,lapply(seq_along(unique(foldid)),function(f) results[[i]]$prediction_errors[[f]][,l]^2))
      })
      # remove standardization of y
      cvraw_lambda = cvraw_lambda*sd_y^2
      # weight cvraw
      if(is.null(prediction_weights)) prediction_weights = rep(1,nobs)
      # group by foldid
      cvob = glmnet::cvcompute(cvraw_lambda,prediction_weights,foldid,rep(length(h2_indexes),max(foldid)))
      cvraw_lambda = cvob$cvraw#-cvraw_lambda_base
      weights = cvob$weights
      N = cvob$N
      cvm_lambda = apply(cvraw_lambda, 2, weighted.mean, w = weights, na.rm = TRUE)
      cvsd_lambda = sqrt(apply(scale(cvraw_lambda, cvm_lambda, FALSE)^2, 2, weighted.mean,
                               w = weights, na.rm = TRUE)/(N - 1))
      list(cvm_lambda=cvm_lambda,cvsd_lambda=cvsd_lambda,cvraw_lambda=cvraw_lambda)
    })
    min_score_index = sapply(cv_results,function(x) which.min(x$cvm_lambda))
    cvms = sapply(cv_results,function(x) x$cvm_lambda)
    cvm = sapply(seq_along(cv_results),function(i) cv_results[[i]]$cvm_lambda[[min_score_index[i]]])
    cvsd = sapply(seq_along(cv_results),function(i) cv_results[[i]]$cvsd_lambda[[min_score_index[i]]])
    min_scores = cvm
  }
  
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
    # min_scores_index_cv = sapply(seq_along(lambda),function(l) {
    #   # go through all h2's, pull out score vectors for each fold: (n(fold) x n(h2) matrix)
    #   lambda_scores_cv = sapply(results,function(x) x$cv_scores[l,])
    #   # for each fold, find the index of the smallest (best) score
    #   apply(lambda_scores_cv,1,which.min)
    # })
    # # pull out prediction_errors for selected models for each fold
    # cvraw = sapply(seq_along(lambda),function(l) {
    #   indexes = min_scores_index_cv[,l]
    #   do.call(c,lapply(seq_along(indexes),function(i) results[[indexes[[i]]]]$prediction_errors[[i]][,l]))^2
    # })
    # # remove standardization of y
    # cvraw = cvraw*sd_y^2
    # # weight cvraw
    # if(is.null(prediction_weights)) prediction_weights = rep(1,nobs)
    # cvm = apply(cvraw, 2, weighted.mean, w = prediction_weights, na.rm = TRUE)
    # cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
    #                   w = prediction_weights, na.rm = TRUE)/(nobs - 1))
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

