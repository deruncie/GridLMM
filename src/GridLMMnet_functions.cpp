#include <math.h>
#include "GridLMM_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;


// [[Rcpp::export()]]
VectorXd Calculate_qt_LASSO(Map<MatrixXd> X, Map<MatrixXd> beta,Map<VectorXd> lambdas){
  // This function calculates tr(X(t(X) %*% X + n*lambda*W^-)^- t(X)), which is the number of effective parameters, part of the calculation of GCV.
  // It uses the fact that W^- is diagonal, with most elements zero, so several low-rank updates can be made to LDLt = t(X) %*% X
  int n = X.rows();
  int p = beta.rows();
  if(X.cols() != p) stop("Wrong dimensions of beta or X");
  int l = lambdas.size();
  if(beta.cols() != l) stop("Wrong length of lambdas");
  MatrixXd nlambda_beta_abs =(1.0/n)*(beta*(lambdas.cwiseInverse().asDiagonal())).cwiseAbs();
  MatrixXd K = X.transpose() * X;
  Eigen::LDLT<MatrixXd,Lower> ldlt_K_base(K);
  VectorXd results = VectorXd::Zero(l);
  for(int i = 0; i < l; i++) {
    Eigen::LDLT<MatrixXd,Lower> ldlt_K = ldlt_K_base;
    for(int j = 0; j < p; j++) {
      if(nlambda_beta_abs.coeffRef(j,i) > 0) {
        VectorXd w = VectorXd::Zero(p);
        w[j] = 1.0/std::sqrt(nlambda_beta_abs.coeffRef(j,i));
        ldlt_K.rankUpdate(w);
      }
    }
    MatrixXd S = X * ldlt_K.solve(X.transpose());
    results[i] = S.diagonal().sum();
  }
  return(results);
}