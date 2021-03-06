% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GridLMMnet.R
\name{GridLMMnet}
\alias{GridLMMnet}
\title{LASSO solutions in a linear mixed model using GridLMM}
\usage{
GridLMMnet(
  formula,
  data,
  X,
  X_ID = "ID",
  weights = NULL,
  centerX = TRUE,
  scaleX = TRUE,
  relmat = NULL,
  normalize_relmat = TRUE,
  h2_step = 0.1,
  h2_start = NULL,
  alpha = 1,
  lambdaType = "s2e",
  scoreType = "LL",
  nlambda = 100,
  lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04),
  lambda = NULL,
  penalty.factor = NULL,
  nfolds = NULL,
  foldid = NULL,
  RE_setup = NULL,
  V_setup = NULL,
  save_V_folder = NULL,
  diagonalize = T,
  mc.cores = parallel::detectCores(),
  clusterType = "mclapply",
  verbose = T,
  ...
)
}
\arguments{
\item{formula}{A two-sided linear formula as used in \code{\link[lme4]{lmer}} describing the fixed-effects
and random-effects of the model on the RHS and the response on the LHS. Note: correlated random-effects
are not implemented, so using one or two vertical bars (\code{|}) or one is identical. At least one random effect is needed.
Unlike \code{lmer}, random effects can have as many as there are observations.}

\item{data}{A data frame containing the variables named in \code{formula}.}

\item{X}{Variables in model that well be penalized with the elastic net penalty. Covariates specified in \code{formula} are not penalized.}

\item{X_ID}{Column of \code{data} that identifies the row of \code{X} that corresponding to each observation. 
It is possible that multiple observations reference the same row of \code{X}.}

\item{weights}{An optional vector of observation-specific weights.}

\item{centerX}{TRUE/FALSE for each. Applied to the \code{X} matrix before using \code{X} to form any GRMs.}

\item{scaleX}{TRUE/FALSE for each. Applied to the \code{X} matrix before using \code{X} to form any GRMs.}

\item{relmat}{Either:
1) A list of matrices that are proportional to the (within) covariance structures of the group level effects. 
2) A list of lists with elements (\code{K}, \code{p}) with a covariance matrix and an integer listing the number of markers 
    used to estimate the covariance matrix. This is used for appropriate downdating of \code{V} to remove proximal markers for each test.
The names of the matrices / list elements should correspond to the columns in \code{data} that are used as grouping factors. All levels
of the grouping factor should appear as rownames of the corresponding matrix.}

\item{normalize_relmat}{should ZKZt matrices be normalized so that mean(diag) == 1? Default (true)}

\item{h2_step}{Step size of the grid}

\item{h2_start}{Optional. Matrix with each row a vector of \eqn{h^2} parameters defining starting values for the grid. 
Typically ML/REML solutions for the null model. 
If null, will be calculated using \link{GridLMM_ML}.}

\item{alpha}{The elasticnet mixing parameter, with \eqn{0\le\alpha\le 1}.
The penalty is defined as
\deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} is the
lasso penalty, and \code{alpha=0} the ridge penalty.}

\item{nlambda}{The number of \code{lambda} values - default is 100.}

\item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of
\code{lambda.max}, the (data derived) entry value (i.e. the smallest value
for which all coefficients are zero). The default depends on the sample size
\code{nobs} relative to the number of variables \code{nvars}. If \code{nobs
> nvars}, the default is \code{0.0001}, close to zero.  If \code{nobs <
nvars}, the default is \code{0.01}.  A very small value of
\code{lambda.min.ratio} will lead to a saturated fit in the \code{nobs <
nvars} case. This is undefined for \code{"binomial"} and
\code{"multinomial"} models, and \code{glmnet} will exit gracefully when the
percentage deviance explained is almost 1.}

\item{lambda}{A user supplied \code{lambda} sequence. Typical usage is to
have the program compute its own \code{lambda} sequence based on
\code{nlambda} and \code{lambda.min.ratio}. Supplying a value of
\code{lambda} overrides this. WARNING: use with care. Avoid supplying a
single value for \code{lambda} (for predictions after CV use
\code{predict()} instead).  Supply instead a decreasing sequence of
\code{lambda} values. \code{glmnet} relies on its warms starts for speed,
and its often faster to fit a whole path than compute a single fit.}

\item{penalty.factor}{Separate penalty factors can be applied to each
coefficient. This is a number that multiplies \code{lambda} to allow
differential shrinkage. Can be 0 for some variables, which implies no
shrinkage, and that variable is always included in the model. Default is 1
for all variables (and implicitly infinity for variables listed in
\code{exclude}). Note: the penalty factors are internally rescaled to sum to
nvars, and the lambda sequence will reflect this change.}

\item{foldid}{vector of integers that divide the data into a set of non-overlapping folds for cross-validation.}

\item{V_setup}{Optional. A list produced by a GridLMM function containing the pre-processed V decompositions for each grid vertex, 
or the information necessary to create this. Generally saved from a previous run of GridLMM on the same data.}

\item{save_V_folder}{Optional. A character vector giving a folder to save pre-processed V decomposition files for future / repeated use. 
If null, V decompositions are stored in memory}

\item{diagonalize}{If TRUE and the model includes only a single random effect, the "GEMMA" trick will be used to diagonalize V. This is done
by calculating the SVD of K, which can be slow for large samples.}

\item{mc.cores}{Number of processor cores used for parallel evaluations. Note that this uses 'mclapply', so the memory requires grow rapidly with \code{mc.cores}, because
the marker matrix gets duplicated in memory for each core.}

\item{verbose}{Should progress be printed to the screen?}

\item{...}{}
}
\value{
If \code{foldid} and \code{nfold} are null, an object with S3 class "glmnet","*" , where "*" is "elnet". See \code{\link[glmnet]{glmnet}}.
    Otherwise, an object with S3 class "cv.glmnet". See \code{\link[glmnet]{cv.glmnet}}.
}
\description{
Finds LASSO or Elastic Net solutions for a multiple regression problem with correlated errors.
}
\details{
Finds the full LASSO or Elastic Net solution path by running \code{\link[glmnet]{glmnet}} at each grid vertex.
    If \code{foldid} is provided, cross-validation scores will be calculated.
}
