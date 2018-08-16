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