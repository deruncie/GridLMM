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
MatrixXd t1(Map<MatrixXd> R,Map<MatrixXd> Y){
  return(R.triangularView<Upper>().solve(Y));
}

// [[Rcpp::export()]]
MatrixXd invt(Map<MatrixXd> Rinv,Map<MatrixXd> Y){
  return(Rinv.triangularView<Upper>() * Y);
}


// [[Rcpp::export()]]
void temp(SEXP Rinv_,Map<MatrixXd> Y, bool sparseR = false){
  if(Rf_isMatrix(Rinv_) ){
    Rcout << "yes";
  } else{
    Rcout << "no";
  }
}
  
// [[Rcpp::export()]]
MatrixXd invt2(SEXP Rinv_,Map<MatrixXd> Y, bool sparseR = false){
  MatrixXd result;
  if(sparseR) {
    MSpMat Rinv = as<MSpMat>(Rinv_);
    result = Rinv.triangularView<Upper>() * Y;
  } else{
    Map<MatrixXd> Rinv = as<Map<MatrixXd> >(Rinv_);
    result = Rinv.triangularView<Upper>() * Y;
  }
  return(result);
}


// [[Rcpp::export()]]
MatrixXd m1(Map<MatrixXd> X, Map<MatrixXd> Y) {
  MatrixXd result(X.rows(),Y.cols());
  result = X*Y;
  return(result);
}
// [[Rcpp::export()]]
MatrixXd m1b(Map<MatrixXd> X, Map<MatrixXd> Y) {
  MatrixXd result(X.rows(),Y.cols());
  for(int i = 0; i < Y.cols(); i++){
    result.col(i) = X*Y.col(i);
  }
  return(result);
}
// [[Rcpp::export()]]
MatrixXd m1c(Map<MatrixXd> Xt, Map<MatrixXd> Y) {
  MatrixXd result(Xt.cols(),Y.cols());
  for(int i = 0; i < Y.cols(); i++){
    result.col(i) = Xt.transpose()*Y.col(i);
  }
  return(result);
}

// [[Rcpp::export()]]
MatrixXd m2(Map<MatrixXd> X, Map<MatrixXd> Y) {
  int n = X.cols();
  int b = Y.rows()/n;
  int p = Y.cols();
  MatrixXd result(Y.rows(),Y.cols());
  for(int i = 0; i < b; i++) {
    result.block(n*i,0,n,p) = X * Y.block(n*i,0,n,p);
  }
  return(result);
}

// [[Rcpp::export()]]
MatrixXd m3(Map<MatrixXd> X, Map<MatrixXd> Y) {
  int n = X.cols();
  int b = Y.rows()/n;
  int p = Y.cols();
  MatrixXd result(Y.rows(),Y.cols());
  for(int i = 0; i < p; i++) {
    MatrixXd Yi = Map<MatrixXd>(Y.col(i).data(),n,b);
    MatrixXd r = X * Yi;
    // result.col(i) = Map<VectorXd>(r.data(),n*b);
  }
  return(result);
}

// [[Rcpp::export()]]
MatrixXd m4(MatrixXd &X, MatrixXd &Y) {
  return(X*Y);
}
// [[Rcpp::export()]]
MatrixXd m4b(MatrixXd X, MatrixXd Y) {
  return(X*Y);
}
