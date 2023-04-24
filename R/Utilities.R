prepMM = function(formula,data,weights = NULL,other_formulas = NULL,
                  relmat = NULL, normalize_relmat = TRUE, X = NULL, X_ID = 'ID',proximal_markers = NULL,V_setup = NULL,
                  diagonalize = TRUE,svd_K = TRUE,drop0_tol = 1e-10,save_V_folder = NULL,
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
  
  # if(!is.null(X_ID) && !X_ID %in% names(RE_terms$cnms)) warning(sprintf("No covariance given SNPs specified in error formula. To specify, add a term like (1|%s)",X_ID))
  
  if(is.null(V_setup)) {
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
      if(term %in% names(relmat)) {
        if(is.list(relmat[[term]])){
          if(verbose && !is.null(proximal_markers)) print('Note: X should be centered and scaled as it was for calculating RRM')
          K = relmat[[term]]$K
          p = relmat[[term]]$p
          p_test = relmat[[term]]$p_test
        } else{
          K = relmat[[term]]
        }
      } else if(!is.null(X_ID) && term == X_ID && !is.null(X)) {
        if(verbose) print('making RRM matrix')
        p = ncol(X)
        K = tcrossprod(X)/p
        relmat[[term]] = K
      } 
      
      if(term == X_ID && !is.null(proximal_markers) && is.null(p_test)){
        p_test = mean(lengths(proximal_markers))
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
    V_setup = make_V_setup(RE_setup,normalize_relmat,weights,diagonalize,svd_K = TRUE,drop0_tol = 1e-10,save_V_folder,verbose)
  } else{
    RE_setup = V_setup$RE_setup
  }
  
  return(list(lmod = lmod, RE_setup = RE_setup, V_setup = V_setup))
}

make_V_setup = function(RE_setup,
                        normalize_relmat = TRUE,
                        weights = NULL,
                        diagonalize = TRUE,
                        svd_K = TRUE, # should the svd of one of the random effect covariances be calculated?
                        drop0_tol = 1e-10, # tolerance for setting a value to zero in a sparse matrix
                        save_V_folder = NULL, # if NULL, V decompositions will not be saved. Otherwise, a string giving the folder to save to. This folder will be cleared of all existing V_decompositions
                        verbose = T
){
  # either returns a folder containing 1) each chol_V file (containing a chol_V, h2_index, and V_log_det) and a chol_Vi_inv_setup file, 2) the h2s_matrix, and 3) the Qt matrix
  # or, a list containing these three items in memory
  
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
  n_RE = length(RE_setup)
  RE_names = names(RE_setup)
  
  
  # ------------------------------------ #
  # ---- Check weights ----------------- #
  # ------------------------------------ #
  
  n = nrow(RE_setup[[1]]$Z)
  if(is.null(weights)) weights = rep(1,n)
  if(length(weights) != n) stop('weights must have same length as data')
  
  
  # ------------------------------------ #
  # ---- Calculate Qt ------------------ #
  # ------------------------------------ #
  
  if(n_RE == 1 && diagonalize == TRUE && all(weights == 1)){
    # only useful if n_RE == 1 and there are no weights
    # it is possible to do this with weights: premultiply Z by diag(sqrt(weights)) below, then set Q = t(U) %*% diag(sqrt(weights)), then adjust log_LL's by sum(log(weights))/2
    # but this is not currently implemented.
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
      if(inherits(ZKZt,'Matrix') && length(ZKZt@x) > ncol(ZKZt)^2/2) ZKZt = as.matrix(ZKZt)
      if(inherits(ZKZt,'Matrix')) {
        result = svd(ZKZt)
      } else{
        result = svd_c(ZKZt)
      }
      Qt = t(result$u)
    }
    # Qt = as(drop0(as(Qt,'dgCMatrix'),tol = drop0_tol),'dgCMatrix')
    QtZ_matrices = lapply(RE_setup,function(re) Qt %*% re$Z)
  } else{
    Qt = NULL
    QtZ_matrices = lapply(RE_setup,function(re) re$Z)
  }
  QtZ = do.call(cbind,QtZ_matrices[RE_names])
  # QtZ = as(QtZ,'dgCMatrix')
  
  # ------------------------------------ #
  # ---- Calculate ZKZts --------------- #
  # ------------------------------------ #
  ZKZts = list()
  for(re in RE_names) {
    # ZKZts[[re]] = as(forceSymmetric(drop0(QtZ_matrices[[re]] %*% RE_setup[[re]]$K %*% t(QtZ_matrices[[re]]),tol = drop0_tol)),'dgCMatrix')
    # ZKZts[[re]] = forceSymmetric(drop0(QtZ_matrices[[re]] %*% RE_setup[[re]]$K %*% t(QtZ_matrices[[re]]),tol = drop0_tol))
    # ZKZts[[re]] = as.matrix(forceSymmetric(drop0(QtZ_matrices[[re]] %*% RE_setup[[re]]$K %*% t(QtZ_matrices[[re]]),tol = drop0_tol)))
    ZKZts[[re]] = QtZ_matrices[[re]] %*% RE_setup[[re]]$K %*% t(QtZ_matrices[[re]])
    if(normalize_relmat) ZKZts[[re]] = ZKZts[[re]]/mean(diag(ZKZts[[re]]))
  }
  
  
  # ------------------------------------ #
  # ---- Prepare setup ------------------ #
  # ------------------------------------ #
  setup = list(
    RE_setup = RE_setup,
    Qt = Qt,
    ZKZts = ZKZts,
    resid_V = diag(1/weights),
    save_V_folder = save_V_folder
  )
  
  # ------------------------------------ #
  # ---- Prep for downdating ----------- #
  # ------------------------------------ #
  setup = set_p_test(setup)
  
  # ------------------------------------ #
  # ---- Prep folder to save V files --- #
  # ------------------------------------ #
  clean_V_folder(setup)
  
  return(setup)
}

clean_V_folder = function(V_setup) {
  save_V_folder = V_setup$save_V_folder
  RE_names = names(V_setup$RE_setup)
  if(!is.null(save_V_folder) & is.character(save_V_folder)){
    try(dir.create(save_V_folder,showWarnings = FALSE),silent = TRUE)
    # clean up existing files
    existing_files = list.files(path = save_V_folder,pattern = 'chol_V_',full.names=TRUE)
    if(length(existing_files) > 0) file.remove(existing_files)
    # prepare index file
    chol_V_index = setNames(data.frame(matrix(ncol = 1+length(RE_names),nrow = 0)),c('file',RE_names))
    write.csv(chol_V_index,file = sprintf('%s/chol_V_index.csv',save_V_folder),row.names=FALSE)
    saveRDS(V_setup,file = sprintf('%s/chol_V_setup.rds',save_V_folder))
  }
}

#' Set the number of markers used in a GRM
#'
#' For GRMs created as X'X, it can be useful to remove "proximal" markers from the GRM when testing a focal marker.
#' To do this in GridLMM, we need to know \code{p} (how many markers were in the original GRM), and \code{p_test},
#' how many markers are to be \strong{removed} for a typical test
#'
#' @param V_setup A list of setup variables, as returned by a GridLMM function
#' @param p_test A named vector with names corresponding to random effects in the model (as found in \code{V_setup$RE_setup}), 
#'    and values equal to the typical number of markers to be \strong{removed} for a typical test
#' @param p A named vector as for \code{p_test} giving the total number of markers used to create the original GRM
#'
#' @return The modified \code{V_setup} list
#' @export
set_p_test = function(V_setup,p_test = NULL, p = NULL){
  RE_setup = V_setup$RE_setup
  if(!is.null(p)) {
    for(re in names(p)){
      if(!re %in% names(RE_setup)) stop(sprintf('Random effect %s not in V_setup',re))
      RE_setup[[re]]$p = p[[re]]
    }
  }
  if(!is.null(p_test)) {
    for(re in names(p_test)){
      RE_setup[[re]]$p_test = p_test[[re]]
    }
  }
  # in preparation for down-dating, multiply each ZKZt by p/(p-pj), where p is the number of SNPs that went into the RRM, and pj is the (mean) number of SNPs per down-date operation
  downdate_ratios = c()
  n_SNPs_downdated_RRM = c()
  for(re in names(RE_setup)) {
    if(is.null(RE_setup[[re]]$p)) {
      downdate_ratios[re] = 1
      n_SNPs_downdated_RRM[re] = Inf
    } else {
      if(is.null(RE_setup[[re]]$p_test)) RE_setup[[re]]$p_test = 0 # assume p_test == 0 if not provided
      n_SNPs_downdated_RRM[re] = RE_setup[[re]]$p - RE_setup[[re]]$p_test
      downdate_ratios[re] = RE_setup[[re]]$p/n_SNPs_downdated_RRM[re]
    }
  }
  V_setup$RE_setup = RE_setup
  V_setup$n_SNPs_downdated_RRM = n_SNPs_downdated_RRM
  V_setup$downdate_ratios = downdate_ratios
  
  # after re-setting p-test, all previously calculated chol_Vi's need to be cleared.
  clean_V_folder(V_setup) 
  
  return(V_setup)
}


make_chol_V_setup = function(V_setup,h2s){
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
  V = (1-sum(h2s)) * V_setup$resid_V
  for(i in 1:length(h2s)) V = V + h2s[i] * as.matrix(V_setup$ZKZts[[i]]) * downdate_ratios[i]
  # test if V is diagonal or sparse
  non_zero_upper_tri_V = abs(V[upper.tri(V,diag = FALSE)]) > 1e-10
  if(sum(non_zero_upper_tri_V) == 0) {  # V is diagonal
    chol_V = as(as(diag(sqrt(diag(V))),'CsparseMatrix'),'dgCMatrix')
  } else if(sum(non_zero_upper_tri_V) < length(non_zero_upper_tri_V) / 2) { # V is sparse
    V = Matrix(V)
    chol_V = as(as(chol(V),'CsparseMatrix'),'dgCMatrix')
  } else{ # V is dense, use Eigen LLT
    chol_V = chol_c(V)
  }
  if(inherits(chol_V,'dtCMatrix')) chol_V = as(chol_V,'dgCMatrix')
  V_log_det = 2*sum(log(diag(chol_V)))
  
  chol_V_setup = list(chol_V = chol_V,V_log_det = V_log_det, h2s = h2s)
  
  if(!is.null(V_setup$save_V_folder)) {
    # save chol_V_setup as a file, append the file to the chol_V_list csv
    chol_V_file = sprintf('%s/chol_V_%s.rds',V_setup$save_V_folder,paste(h2s,collapse='_'))
    saveRDS(chol_V_setup,file = chol_V_file)
    system(sprintf('echo "%s,%s" >> %s/chol_V_index.csv',chol_V_file,paste(h2s,collapse=','),V_setup$save_V_folder))
    chol_V_setup$file = chol_V_file
  }
  return(chol_V_setup)
}


scale_SNPs = function(X,centerX=TRUE,scaleX=TRUE,fillNAX = FALSE){
  if(!centerX & !scaleX & !fillNAX) {
    if(any(is.na(X))) stop('NAs in X. Set fillNAX = TRUE')
    return(X)
  }
  if(all(range(X,na.rm=T) == c(-1,1))){
    pi = colMeans((X+1)/2,na.rm=T)
    # if(centerX){
    #   X = sweep(X,2,2*pi,'-')
    # }
    if(scaleX){
      X = sweep(X,2,sqrt(2*pi*(1-pi)),'/')
    }
  } else if(all(range(X,na.rm=T) == c(0,1))){
    pi = colMeans(X,na.rm=T)
    if(centerX){
      X = sweep(X,2,pi,'-')
    }
    if(scaleX){
      X = sweep(X,2,sqrt(pi*(1-pi)),'/') # is this pq?
    }
  } else if(all(range(X,na.rm=T) == c(0,2))){
    pi = colMeans(X/2,na.rm=T)
    if(centerX){
      X = sweep(X,2,2*pi,'-')
    }
    if(scaleX){
      X = sweep(X,2,sqrt(2*pi*(1-pi)),'/')
    }
  } else{
    if(centerX){
      X = sweep(X,2,colMeans(X,na.rm=T),'-')
    }
    if(scaleX) {
      X = sweep(X,2,apply(X,2,sd,na.rm=T),'/')
    }
  }
  if(sum(is.na(X))>0){
    if(fillNAX){
      X = apply(X,2,function(x) {
        if(sum(is.na(x)) > 0){
          x[is.na(x)] = mean(x,na.rm=T)
        }
        x
      })
    } else{
      X[,colSums(is.na(X))>1] = 0
    }
  }
  X
}

# makes a grid over a set of h2s.
# very simlar to make_h2s_matrix
# probably could/should combine functions
setup_Grid = function(RE_names,h2_step,h2_start = NULL){
  
  # -------- Form matrix of h2s to test ---------- #
  n_RE = length(RE_names)
  
  if(length(h2_step) < n_RE){
    if(length(h2_step) != 1) stop('Must provide either 1 h2_step parameter, or 1 for each random effect')
    h2_step = rep(h2_step,n_RE)
  }
  if(is.null(names(h2_step))) {
    names(h2_step) = RE_names
  }
  
  if(is.null(h2_start)) h2_start = rep(0,length(RE_names))
  if(length(h2_start) != length(RE_names)) stop("Wrong length of h2_start")
  if(is.null(names(h2_start))) names(h2_start) = RE_names
  
  h2s_matrix = expand.grid(lapply(RE_names,function(re) h2_start[[re]] + seq(-1,2,by = h2_step[[re]])))
  colnames(h2s_matrix) = RE_names
  h2s_matrix = h2s_matrix[rowSums(h2s_matrix<0) == 0,,drop=FALSE]
  h2s_matrix = t(h2s_matrix[rowSums(h2s_matrix) < 1,,drop=FALSE])
  colnames(h2s_matrix) = NULL
  
  return(h2s_matrix)
}

calculate_Grid = function(V_setup,h2s_matrix,verbose = TRUE){
  
  if(verbose) {
    sprintf('Generating V decompositions for %d grid cells', ncol(h2s_matrix))
    pb = txtProgressBar(min=0,max = ncol(h2s_matrix),style=3)
  }
  
  V_setup$chol_V_list = foreach(h2s = iter(t(h2s_matrix),by='row')) %dopar% {
    h2s = h2s[1,]
    chol_V_setup = make_chol_V_setup(V_setup,h2s)
    if(!is.null(V_setup$save_V_folder)) chol_V_setup = chol_V_setup$file  # only store file name if save_V_folder is provided
    if(verbose && exists('pb')) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    return(chol_V_setup)
  }
  if(verbose && exists('pb')) close(pb)
  
  return(V_setup)
}



# makes a grid over a set of h2s
# checks the grid for valid h2s_vectors (rows)
# RE_names is the names of each element of the h2s vector
# steps is the grid values for each h2 element
make_h2s_matrix = function(RE_names,steps) {
  if(!is.list(steps)) {
    steps = lapply(RE_names,function(x) steps)
    names(steps) = RE_names
  }
  if(length(steps) != length(RE_names)) stop('wrong length of grid points. If a list, should have the length of RE_names')
  if(is.null(names(steps))) names(steps) = RE_names
  
  h2s_matrix = expand.grid(steps[RE_names])
  colnames(h2s_matrix) = make.names(RE_names)
  h2s_matrix = as.matrix(h2s_matrix[rowSums(h2s_matrix) < 1-1e-10,,drop=FALSE])
  rownames(h2s_matrix) = NULL
  h2s_matrix
}


# collects all unique h2s values from a results data.frame
# h2s columns are named: variable.ML or variable.REML
# returns a matrix with rows as h2s vectors
get_current_h2s = function(results,RE_names,ML,REML) {
  rownames(results) = NULL
  current_h2s = c()
  if(ML) {
    current_h2s_ML = unique(as.matrix(results[,paste(make.names(RE_names),'ML',sep='.'),drop=FALSE]))
    colnames(current_h2s_ML) = sub('.ML','',colnames(current_h2s_ML),fixed=T)
    current_h2s = current_h2s_ML
  }
  if(REML){
    current_h2s_REML = unique(as.matrix(results[,paste(make.names(RE_names),'REML',sep='.'),drop=FALSE]))
    colnames(current_h2s_REML) = sub('.REML','',colnames(current_h2s_REML),fixed=T)
    current_h2s = unique(rbind(current_h2s,current_h2s_REML))
  }
  current_h2s = na.omit(current_h2s) # some tests fail
  return(current_h2s)
}

# makes balls around each row of a h2s_matrix
# checks the balls for valid vectors
get_h2s_ball = function(current_h2s,h2_step){
  RE_names = colnames(current_h2s)
  h2s_ball = make_h2s_matrix(RE_names,c(-1,0,1)*h2_step)
  
  if(nrow(h2s_ball) == 0) {
    print(current_h2s)
    print(h2_step)
    print(h2s_ball)
  }
  
  h2s_to_test = c()
  for(i in 1:nrow(current_h2s)){
    h2s = current_h2s[i,]
    new_h2s = sweep(h2s_ball,2,unlist(h2s),'+')
    h2s_to_test = rbind(h2s_to_test,new_h2s)
  }
  h2s_to_test = unique(h2s_to_test)
  h2s_to_test = h2s_to_test[apply(h2s_to_test,1,min) >= 0,,drop=FALSE]
  h2s_to_test = h2s_to_test[apply(h2s_to_test,1,max) < 1,,drop=FALSE]
  h2s_to_test = h2s_to_test[rowSums(h2s_to_test) < 1-1e-10,,drop=FALSE]
  rownames(h2s_to_test) = NULL
  return(h2s_to_test)
}


# finds current h2s, makes balls around them, compares them to tested_h2s, returns matrix
# the challenge is that for some reason get_h2s_ball gets off by a very small amount, so h2s vectors aren't identical
get_h2s_to_test = function(current_h2s,tested_h2s,h2_step,ML,REML){ 
  all_h2s_to_test = get_h2s_ball(current_h2s,h2_step)
  h2s_to_test = tested_h2s
  for(i in 1:nrow(all_h2s_to_test)){
    n = nrow(h2s_to_test)
    h2s = all_h2s_to_test[i,]
    dists = sqrt(rowSums(sweep(h2s_to_test,2,h2s,'-')^2))
    if(min(dists[1:n]) >= 1e-10) {
      h2s_to_test = rbind(h2s_to_test,h2s)
    } 
  }
  h2s_to_test = h2s_to_test[-c(1:nrow(tested_h2s)),,drop=FALSE]
  rownames(h2s_to_test) = NULL
  return(h2s_to_test)
}

calc_ML = function(SSs,n) {
  return(n/2*log(n/(2*pi)) - n/2 - SSs$V_log_dets/2 - n/2*log(t(SSs$RSSs)))
}

calc_REML = function(SSs,n,b,m){
  return(SSs$ML + 1/2*(b*log(2*pi*t(SSs$RSSs)/(n-b)) + SSs$log_det_X - SSs$V_star_inv_log_det)) # This is the formula from Kang 2008. It's also the same in FaST-LMM.
  # df_E = n - b # This is the formula from Gemma 2012
  # return(df_E/2*log(df_E/(2*pi)) - df_E/2 + SSs$log_det_X/2 - SSs$V_log_dets/2 - SSs$V_star_inv_log_det/2 - df_E/2*t(log(SSs$RSSs)))
}

get_LL = function(SSs,X_cov,X_list,active_X_list,n,m,ML,REML,BF){
  SSs$ML = calc_ML(SSs,n)
  SSs$s2 = SSs$RSSs / n
  if(REML) {
    if(is.null(SSs$log_det_X)) SSs$log_det_X = log_det_of_XtX(X_cov,X_list,active_X_list)
    b = ncol(X_cov) + length(X_list)
    if(length(X_list) == 0) b = ncol(X_cov)
    # The following is not needed. Should use eq 13 from Kang et al 2008. Note: if M is the identity with the first k rows removed, (M(X'VinvX)M')' is just crossprod(V_star_L[-c(1:k),-c(1:k)])
    # somehow need to figure out how to specify M
    SSs$F_hats = F_hats(SSs$beta_hats,SSs$RSSs,SSs$V_star_L,n,b,m)  
    SSs$REML = calc_REML(SSs,n,b,m)
    SSs$s2 = SSs$RSSs / (n-b)
  }
  if(BF){
    a_star = n/2
    b_star = SSs$RSSs/2
    SSs$log_posterior_factor = - SSs$V_star_inv_log_det/2 - SSs$V_log_dets/2 - a_star*log(b_star)
  }
  return(SSs)
}


calc_LL = function(Y,X_cov,X_list,h2s,chol_Vi,inv_prior_X,
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
    SSs <- GridLMM_SS_matrix(Y,chol_Vi,X_cov,X_list,active_X_list,inv_prior_X)
  } else{
    downdate_weights = h2s/unlist(n_SNPs_downdated_RRM)
    # the length of downdate_Xi determines how many K's get downdated
    # this means that h2s has to be ordered such that the downdated K's are first and in the correct order
    downdate_weights = downdate_weights[1:length(downdate_Xs[[active_X_list[[1]]]]$downdate_Xi)] 
    if(length(downdate_weights) != length(downdate_Xs[[active_X_list[[1]]]]$downdate_Xi)) stop("Wrong length of downdate weights")
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
                         s2 = c(t(log_LL$s2)),
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

calc_LL_parallel = function(Y,X_cov,X_list,h2s,chol_V,inv_prior_X,
                             downdate_Xs = NULL,n_SNPs_downdated_RRM = NULL,REML = TRUE,BF = FALSE,
                             mc.cores = 1,active_X_list = NULL) {
  
  if(mc.cores == 1 || length(X_list) == 0) {
    return(calc_LL(Y,X_cov,X_list,h2s,chol_V,inv_prior_X,downdate_Xs,n_SNPs_downdated_RRM,REML,BF,active_X_list))
  }
  # make X_list_active, divide into chunks
  X_list_active = NULL
  downdate_Xs_active = NULL
  total_tests = length(active_X_list)
  if(!is.null(X_list)) {
    X_list_active = lapply(X_list,function(x) x[,active_X_list,drop=FALSE])
  }
  if(!is.null(downdate_Xs)) {
    downdate_Xs_active = downdate_Xs[active_X_list]
  }
  
  # X_list_names = colnames(X_list[[1]])
  registerDoParallel(mc.cores)
  chunkSize = total_tests/mc.cores
  chunks = 1:total_tests
  chunks = split(chunks, ceiling(seq_along(chunks)/chunkSize))
  X_list_sets = lapply(chunks,function(i) lapply(X_list_active,function(Xi) Xi[,i,drop=FALSE]))
  if(!is.null(downdate_Xs_active)) {
    downdate_Xs_sets = lapply(chunks,function(i) downdate_Xs_active[i])
    results = foreach(X_list_i = iter(X_list_sets),downdate_Xs_i = iter(downdate_Xs_sets),.combine = 'rbind') %dopar% {
      calc_LL(Y,X_cov,X_list_i,h2s,chol_V,inv_prior_X,downdate_Xs_i,n_SNPs_downdated_RRM,REML,BF)
    }
  } else{
    results = foreach(X_list_i = iter(X_list_sets),.combine = 'rbind') %dopar% {
      calc_LL(Y,X_cov,X_list_i,h2s,chol_V,inv_prior_X,downdate_Xs,n_SNPs_downdated_RRM,REML,BF)
    }
  }
  return(results)
}


# takes a list of data.frames of results for the same tests and chooses the best value by ML (and REML), returning a single data.frame of results
compile_results = function(results_list){
  results = results_list[[1]]
  if(length(results_list) == 1) return(results_list[[1]])
  
  results1_ID = paste(results$Trait,results$X_ID,sep='::')
  results_list_ID = lapply(results_list,function(x) {
    id = paste(x$Trait,x$X_ID,sep='::')
    match(results1_ID,id)
  })
  
  ML_h2_columns = colnames(results)[grep(".ML",colnames(results),fixed=T)]
  beta_columns = colnames(results)[grep("beta.",colnames(results),fixed=T)]
  REML_h2_columns = colnames(results)[grep(".REML",colnames(results),fixed=T)]
  F_columns = colnames(results)[grep("F.",colnames(results),fixed=T)]
  posterior_column = colnames(results)[grep("log_posterior_factor",colnames(results),fixed=T)]
  tau_column = colnames(results)[grep("Tau",colnames(results),fixed=T)]
  s2_column = colnames(results)[grep("s2",colnames(results),fixed=T)]
  
  ML = FALSE
  if(length(ML_h2_columns) > 0) ML = TRUE
  REML = FALSE
  if(length(REML_h2_columns) > 0) REML = TRUE
  BF = FALSE
  if(length(posterior_column) > 0) BF = TRUE
  
  if(ML) {
    MLs <- foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i$ML_logLik[index] 
    # max_ML_index = apply(MLs,1,function(x) which.max(x)[1])
    MLs[is.na(MLs)] = -Inf
    max_ML_index = max.col(MLs,'first')
    results$ML_index = max_ML_index
    
    ML_index = 1:nrow(MLs) + (max_ML_index-1)*nrow(MLs) # this pulls out the ML value from each row from the column given by max_ML_index
    results$ML_logLik = MLs[ML_index] 
    
    for(ML_h2_column in ML_h2_columns) {
      ML_h2s = foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i[index,ML_h2_column]
      results[[ML_h2_column]] = ML_h2s[ML_index]
    }
  }
  
  if(!REML) {
    for(beta_column in beta_columns) {
      betas = foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i[index,beta_column]
      results[[beta_column]] = betas[ML_index]
    }
  } else{
    REMLs <- foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i$REML_logLik[index] 
    # max_REML_index = apply(REMLs,1,function(x) which.max(x)[1])
    REMLs[is.na(REMLs)] = -Inf
    max_REML_index = max.col(REMLs,'first')
    results$REML_index = max_REML_index
    
    REML_index = 1:nrow(REMLs) + (max_REML_index-1)*nrow(REMLs) # this pulls out the REML value from each row from the column given by max_ML_index
    results$REML_logLik = REMLs[REML_index] 
    
    for(REML_h2_column in REML_h2_columns) {
      REML_h2s = foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i[index,REML_h2_column]
      results[[REML_h2_column]] = REML_h2s[REML_index]
    }
    
    for(beta_column in beta_columns) {
      betas = foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i[index,beta_column]
      results[[beta_column]] = betas[REML_index]
    }
    
    for(F_column in F_columns) {
      F_hats = betas = foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i[index,F_column]
      results[[F_column]] = F_hats[REML_index]
    }
    
    if(length(tau_column)>0) {
      taus = foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i[index,tau_column]
      results[[tau_column]] = taus[REML_index]
    }
    
    if(length(s2_column)>0) {
      s2s = foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i[index,s2_column]
      results[[s2_column]] = s2s[REML_index]
    }
    
  }
  
  
  if(BF) {
    log_posterior_factors = foreach(results_i = results_list,index = results_list_ID,.combine = 'cbind') %do% results_i[index,posterior_column]
    # log_posterior_factors = foreach(x = iter(log_posterior_factors,by='row'),.combine = 'rbind') %do% sort(x,decreasing=T,na.last=T)
    if(is.null(dim(log_posterior_factors))) log_posterior_factors = matrix(log_posterior_factors,nr=1)
    # add posterior factors without overflow - keep in log space
    max_log_posterior_factors = apply(log_posterior_factors,1,max,na.rm=T)
    log_total_posterior_factors = max_log_posterior_factors + log(rowSums(exp(log_posterior_factors - max_log_posterior_factors),na.rm=T))
    results[[posterior_column]] = log_total_posterior_factors
    #
    # calculate the amount of posterior mass contributed by the new grid locations
    # need to add the prior here!
    posteriors = exp(log_posterior_factors - log_total_posterior_factors)
    sum_new_posterior = rowSums(posteriors[,-1,drop=FALSE],na.rm=T)
    results$BF_new_posterior_mass = sum_new_posterior
  }
  
  return(results)
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

stop_cluster = function(cl){
  try(parallel::stopCluster(cl))
  foreach::registerDoSEQ()
}

#' Prepare a Kinship matrix for use with LDAK
#'
#' Takes an R matrix and saves it to disk in a format that LDAK can use
#' 
#' @details The matrix is first saved as a text file (without row/column names), and then 
#'     LDAK is called by \code{system} with option "--convert-raw" to convert it into the format
#'     used by LDAK. Temporary files are removed. The .grn.id file contains the rownames of \code{K} repeated in two columns.
#' 
#' @param K A matrix. Must have rownames.
#' @param kinfile The base of the filename fo the kinship files to be created, including path.
#' @param LDAK_program The path to the LDAK program
#'
#' @return A character with the base file name of the converted matrix.
#' @export
#'
prep_LDAK_Kinship = function(K,kinfile,LDAK_program = 'LDAK/ldak5') {
  # data.table::fwrite(data.frame(as.matrix(K)),file = 'temp.grm.raw',row.names=F,col.names=F,quote=F)
  write.table(as.matrix(K),file = 'temp.grm.raw',row.names=F,col.names=F,quote=F)
  write.table(cbind(rownames(K),rownames(K)),file = 'temp.grm.id',row.names=F,col.names=F,quote=F)
  system(sprintf('./%s --convert-raw %s --grm temp',LDAK_program,kinfile))
  system(sprintf('rm temp.grm.*',kinfile))
  return(kinfile)
}

combine_LDAK_Kinships = function(K_list, file = 'K.list') {
  kinship_list = c('trait_A_converted','trait_D_converted','trait_E_converted','cage_K_converted','pop_K_converted')
  write.table(K_list,file = file,row.names=F,col.names=F,quote=F)
  return(file)
}

prep_h2_LDAK = function(y,X_cov,K_list,LDAK_program, maxiter = 1000){
  ID = data.table::fread(paste0(K_list[1],'.grm.id'),data.table = F,h=F)[,1:2]
  write.table(cbind(ID,y),file = 'phen.txt',row.names=F,col.names=F,quote=F)
  if(all(X_cov[,1] == 1)) X_cov = X_cov[,-1,drop=FALSE]
  write.table(cbind(ID,X_cov),file = 'cov.txt',row.names=F,col.names=F,quote=F)
  combine_LDAK_Kinships(K_list, file = 'K.list')
  if(ncol(X_cov) > 0) {
    return(sprintf('./%s --reml h2_LDAK --pheno phen.txt --covar cov.txt --mgrm K.list --kinship-details NO --constrain YES --reml-iter %d',LDAK_program,maxiter))
  } else{
    return(sprintf('./%s --reml h2_LDAK --pheno phen.txt --mgrm K.list --kinship-details NO --constrain YES --reml-iter %d',LDAK_program,maxiter))
  }
}
get_LDAK_result = function(K_list,weights = rep(1,length(K_list))){
  ldak_results = fread('h2_LDAK.vars',data.table = F)
  ldak_results$Variance[-nrow(ldak_results)] = ldak_results$Variance[-nrow(ldak_results)]/weights
  h2_LDAK = matrix(ldak_results$Variance[-nrow(ldak_results)]/sum(ldak_results$Variance),nr=1)
  colnames(h2_LDAK) = names(K_list)
  rownames(h2_LDAK) = NULL
  h2_LDAK
}
get_h2_LDAK = function(y,X_cov,K_list,LDAK_program,weights = rep(1,length(K_list)), maxiter = 1000){
  LDAK_call = prep_h2_LDAK(y,X_cov,K_list,LDAK_program,maxiter)
  system(LDAK_call)
  get_LDAK_result(K_list,weights)
}

time_LDAK = function(y,X_cov,K_list,LDAK_program){
  times = c()
  for(iters in c(0,1,2,4,6,10,1000)){
    LDAK_call = prep_h2_LDAK(y,X_cov,K_list,LDAK_program,iters)
    # LDAK_call = paste(LDAK_call,'--tolerance 0.000000001')
    time_i = system.time({system(LDAK_call)})
    res = fread('h2_LDAK.progress')
    n_iter = nrow(res)
    times = rbind(times,data.frame(n_iter = n_iter,time = time_i[3]))
  }
  times
}


qq_p_plot = function(ps,...){
  if(is.list(ps)) {
    ps_list = lapply(ps,function(x) {
      list(ps = -sort(log10(x)), us = -log10(qunif(ppoints(length(x)))))
    })
    names(ps_list) = names(ps)
    max_p = max(sapply(ps_list,function(x) x$ps[1]))
    max_u = max(sapply(ps_list,function(x) x$us[1]))
    ylim = c(0,max_p)
    xlim = c(0,max_u)
    plot(NA,NA,xlim = xlim,ylim = ylim,xlab = 'Expected -log10(p)',ylab = 'Observed -log10(p)',...)
    abline(0,1)
    for(i in 1:length(ps_list)){
      points(ps_list[[i]]$us,ps_list[[i]]$ps,pch=19,cex=.5,col=i)
    }
    if(!is.null(names(ps))[1]) {
      legend('topleft',legend = names(ps_list),col = 1:length(ps_list),pch=19)
    }
  } else{
    ps = sort(ps)
    u = sort(runif(length(ps)))
    plot(-log10(u),-log10(ps),xlab = 'Expected',ylab = 'Observed',...)
    abline(0,1)
  }
}

