#' Detect number of available cores on system
#'
#' \code{my_detectCores()} detects the numer of cores available for calls to \code{\link{mclapply}}
#' or \code{\link{registerDoParallel}}.
#' 
#' The command is specific for clusters running slurm. 
#' Accesses the \code{SLURM_CPUS_PER_TASK} environmental variable. If this
#' variable doesn't exist it uses \code{\link[parallel]{detectCores}}.
#' 
#' @return ncores integer of number of avaialble cores
#'
#' @export
#'
#' @examples
#' 
#' my_detectCores()
my_detectCores = function() {
  ncores = suppressWarnings(as.numeric(system('printenv SLURM_CPUS_PER_TASK',intern=T)))
  if(length(ncores) == 0 || is.na(ncores)) ncores = parallel::detectCores()
  ncores
}
