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
MatrixXd chol_update(MatrixXd L, MatrixXd X, int sign) {
  int n = L.rows();
  Map<MatrixXd> L_new(L.data(),n,n);
  for(int i = 0; i < X.cols(); i++){
    for(int k = 0; k < n; k++){
      double L_kk = L_new.coeffRef(k,k);
      double x_k = X.coeffRef(k,i);
      double r = sqrt(L_kk*L_kk + sign * x_k*x_k);
      double c = r / L_kk;
      double s = x_k / L_kk;
      L_new.coeffRef(k,k) = r;
      if(k < (n-1)) {
        L_new.block(k+1,k,n-k-1,1) = (L.block(k+1,k,n-k-1,1) + sign * s * X.block(k+1,i,n-k-1,1))/c;
        X.block(k+1,i,n-k-1,1) = c*X.block(k+1,i,n-k-1,1) - s*L_new.block(k+1,k,n-k-1,1);
      }
    }
  }
  return(L_new);
}


// [[Rcpp::export()]]
MatrixXd chol_dropRows(Map<MatrixXd> L,int start_row, int num_rows){
  // go to zero-index
  start_row -= 1;
  
  int n = L.rows();
  int n_new = n - num_rows;
  int num_bottom_rows = n - (start_row + num_rows);
  
  MatrixXd Lnew = MatrixXd::Zero(n_new,n_new);
  
  if(start_row > 0){
    Lnew.topLeftCorner(start_row,start_row) = L.topLeftCorner(start_row,start_row);
  }
  if(num_bottom_rows > 0) {
    Lnew.block(start_row,0,num_bottom_rows,start_row) = L.block(start_row+num_rows,0,num_bottom_rows,start_row);
    
    MatrixXd S32 = L.block(start_row+num_rows,start_row,num_bottom_rows,num_rows);
    MatrixXd S33 = L.bottomRightCorner(num_bottom_rows,num_bottom_rows);
    Lnew.bottomRightCorner(num_bottom_rows,num_bottom_rows) = chol_update(S33,S32,1);
  }
  
  return(Lnew);
}
