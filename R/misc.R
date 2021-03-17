#' Detect number of available cores on system
#'
#' \code{my_detectCores()} detects the number of cores available for calls to \code{\link{mclapply}}
#' or \code{\link{registerDoParallel}}.
#' 
#' This uses RcppParallel::defaultNumThreads() which seems to generally work for me.
#' 
#' @return ncores integer of number of available cores
#'
#' @export
#'
#' @examples
#' 
#' my_detectCores()
my_detectCores = function() {
  ncores = RcppParallel::defaultNumThreads()
  if(length(ncores) == 0 || is.na(ncores)) ncores = parallel::detectCores()
  ncores
}
