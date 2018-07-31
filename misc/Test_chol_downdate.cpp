#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;


using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::ArrayXd;                  // variable size vector, double precision
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;

using namespace Eigen;


// [[Rcpp::export()]]
MatrixXd multiply(Map<MatrixXd> X, Map<MatrixXd> Y) {
  return(X*Y);
}

// [[Rcpp::export()]]
MatrixXd crossprod_cholR(Map<MatrixXd> chol_R, Map<MatrixXd> X){
  // t(chol_R) %*% chol_R = K
  // returns t(chol_R) %*% X
  return(chol_R.transpose().triangularView<Lower>() * X);
}