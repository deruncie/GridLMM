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

void chol_update_R_inplace2(MatrixXd &R, MatrixXd X, VectorXd weights) {
  int n = R.rows();
  if(X.rows() != n) stop("Wrong dimension of X for downdating R");
  if(weights.size() != X.cols()) stop("wrong length of weights for downdating R");
  for(int i = 0; i < X.cols(); i++){
    VectorXd Xi = X.col(i);
    double weights_i = weights[i];
    for(int k = 0; k < n; k++){
      double R_kk = R.coeffRef(k,k);
      double x_k = Xi.coeffRef(k);
      double r = sqrt(R_kk*R_kk + weights_i * x_k*x_k);
      double c = r / R_kk;
      double s = x_k / R_kk;
      R.coeffRef(k,k) = r;
      if(k < (n-1)) {
        for(int j=k+1; j<n; j++) {
          double R_kj = R.coeffRef(k,j);
          double X_ij = Xi.coeff(j);
          R_kj += weights_i * s * X_ij;
          R_kj /= c;
          X_ij *= c;
          X_ij -= s*R_kj;
          R.coeffRef(k,j) = R_kj;
          Xi.coeffRef(j) = X_ij;
        }
        // R.block(k,k+1,1,n-k-1) = (R.block(k,k+1,1,n-k-1) + weights_i * s * Xi.tail(n-k-1).transpose())/c;
        // Xi.tail(n-k-1) = c*Xi.tail(n-k-1) - s*R.block(k,k+1,1,n-k-1).transpose();
        // R.block(k,k+1,1,n-k-1) += weights_i * s * Xi.tail(n-k-1).transpose();
        // R.block(k,k+1,1,n-k-1) /= c;
        // Xi.tail(n-k-1) *= c;
        // Xi.tail(n-k-1) -= s*R.block(k,k+1,1,n-k-1).transpose();
        }
    }
  }
}


// [[Rcpp::export()]]
MatrixXd chol_update_L2(MatrixXd L, MatrixXd X, VectorXd weights) {
  MatrixXd R = L.transpose();
  chol_update_R_inplace2(R,X,weights);
  return(R.transpose());
}

// [[Rcpp::export()]]
MatrixXd chol_update_R2(MatrixXd R, MatrixXd X, VectorXd weights) {
  chol_update_R_inplace2(R,X,weights);
  return(R);
}

void chol_update_L_inplace2(MatrixXd &L, MatrixXd X, VectorXd weights) {
  int n = L.rows();
  if(X.rows() != n) stop("Wrong dimension of X for downdating L");
  if(weights.size() != X.cols()) stop("wrong length of weights for downdating L");
  for(int i = 0; i < X.cols(); i++){
    VectorXd Xi = X.col(i);
    double weights_i = weights[i];
    for(int k = 0; k < n; k++){
      double L_kk = L.coeffRef(k,k);
      double x_k = Xi.coeffRef(k);
      double r = sqrt(L_kk*L_kk + weights_i * x_k*x_k);
      double c = r / L_kk;
      double s = x_k / L_kk;
      L.coeffRef(k,k) = r;
      if(k < (n-1)) {
        for(int j=k+1; j<n; j++) {
          double L_jk = L.coeffRef(j,k);
          double X_ij = Xi.coeff(j);
          L_jk += weights_i * s * X_ij;
          L_jk /= c;
          X_ij *= c;
          X_ij -= s*L_jk;
          L.coeffRef(j,k) = L_jk;
          Xi.coeffRef(j) = X_ij;
        }
        // L.block(k+1,k,n-k-1,1) = (L.block(k+1,k,n-k-1,1) + weights_i * s * Xi.tail(n-k-1))/c;
        // Xi.tail(n-k-1) = c*Xi.tail(n-k-1) - s*L.block(k+1,k,n-k-1,1);
        // R.block(k,k+1,1,n-k-1) += weights_i * s * Xi.tail(n-k-1).transpose();
        // R.block(k,k+1,1,n-k-1) /= c;
        // Xi.tail(n-k-1) *= c;
        // Xi.tail(n-k-1) -= s*R.block(k,k+1,1,n-k-1).transpose();
      }
    }
  }
}

// [[Rcpp::export()]]
MatrixXd chol_update_L3(MatrixXd L, MatrixXd X, VectorXd weights) {
  chol_update_L_inplace2(L,X,weights);
  return(L);
}

// [[Rcpp::export()]]
MatrixXd chol_solve1(Map<MatrixXd> chol_R, Map<MatrixXd> Y){
  return(chol_R.transpose().triangularView<Lower>().solve(Y));
}
// [[Rcpp::export()]]
MatrixXd chol_solve2(Map<MatrixXd> chol_L, Map<MatrixXd> Y){
  return(chol_L.triangularView<Lower>().solve(Y));
}

// [[Rcpp::export()]]
MatrixXd chol_u1(Map<MatrixXd> A, VectorXd X, int times){
  LLT<MatrixXd> lA(A);
  MatrixXd L = lA.matrixL();
  // VectorXd weights = VectorXd::Constant(1,1.0);
  // for(int i = 0; i < times; i++){
  //   chol_update_L_inplace2(L, X, weights);
  // }
  MatrixXd XX = X.replicate(1,times);
  VectorXd weights = VectorXd::Constant(times,1.0);
  chol_update_L_inplace2(L,XX,weights);
  return(L);
}

// [[Rcpp::export()]]
MatrixXd chol_u2(Map<MatrixXd> A, VectorXd X, int times){
  LLT<MatrixXd> lA(A);
  VectorXd weights = VectorXd::Constant(1,1.0);
  for(int i = 0; i < times; i++){
    lA.rankUpdate(X,1.0);
  }
  MatrixXd L = lA.matrixL();
  return(L);
}
// 
// Rcpp::List chol_cT(Map<MatrixXd> A){
//   LLT<MatrixXd> lA(A);
//   return(Rcpp::List::create(lA));
// }
// // [[Rcpp::export()]]
// MatrixXd chol_cL(Rcpp::List L){
//   LLT<MatrixXd> lA = as<LLT<MatrixXd>>(L[0]);
//   return(lA.matrixL());
// }


// void a1(MatrixXd &X){
//   X.coeffRef(0,0) += 1;
//   Rcout << X.coeff(0,0) << std::endl;
// }
// 
// // [[Rcpp::export()]]
// void a0(MatrixXd X){
//   Rcout << X.coeff(0,0) << std::endl;
//   a1(X);
//   Rcout << X.coeff(0,0) << std::endl;
//   MatrixXd Y(X);
//   a1(Y);
//   Rcout << X.coeff(0,0) << std::endl;
// }
// 
// struct x{
//   double A=0;
// };
// 
// x a3a(double A){
//   x X;
//   X.A = A;
//   return(X);
// }
// void a3(x &X){
//   // X.A = 3;
//   X = a3a(5);
// }
// 
// // [[Rcpp::export()]]
// void a4(){
//   x X;
//   Rcout << X.A << std::endl;
//   a3(X);
//   Rcout << X.A << std::endl;
// }

// // [[Rcpp::export()]]
// Rcpp::List svd_c(Map<MatrixXd> X) {
//   Eigen::BDCSVD<MatrixXd> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
//   // bcdsolve.computeU();
//   // MatrixXd U = bcdsolve.matrixU();
//   return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
//                      Named("d") = bcdsolve.singularValues()));
// }
// 
// 
// // [[Rcpp::export()]]
// Rcpp::List svd_c3(SEXP X) {
//   if(Rf_isMatrix(X)) {
//     Eigen::BDCSVD<MatrixXd> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
//     // bcdsolve.computeU();
//     // MatrixXd U = bcdsolve.matrixU();
//     return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
//                               Named("d") = bcdsolve.singularValues()));
//   } else{
//     Eigen::BDCSVD<MSpMat> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
//     // bcdsolve.computeU();
//     // MatrixXd U = bcdsolve.matrixU();
//     return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
//                               Named("d") = bcdsolve.singularValues()));
//   }
// }
// 
// // [[Rcpp::export()]]
// Rcpp::List svd_c2(MatrixXd X) {
//   Eigen::BDCSVD<MatrixXd> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
//   // bcdsolve.computeU();
//   // MatrixXd U = bcdsolve.matrixU();
//   return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
//                             Named("d") = bcdsolve.singularValues()));
// }
// 
// 
// // [[Rcpp::export()]]
// void m1(Map<MatrixXd> X){
//   X *= 3;
// }
// 
// // [[Rcpp::export()]]
// void m1b(Map<MatrixXd> X){
//   MatrixXd X2(X);
//   X2 *= 3;
// }
// 
// 
// // [[Rcpp::export()]]
// void test_list(Rcpp::List X_){
//   std::vector<Map<MatrixXd> > X;
//   for(int i = 0; i < X_.length(); i++){
//     X.push_back(as<Map<MatrixXd> >(X_[i]));
//   }
//   X[0] *= 3;
// }
// // [[Rcpp::export()]]
// void test_list2(Rcpp::List X_){
//   std::vector<MatrixXd> X;
//   for(int i = 0; i < X_.length(); i++){
//     X.push_back(MatrixXd());
//     MatrixXd Xi = as<MatrixXd>(X_[i]);
//     X.back().swap(as<Map<MatrixXd> >(X_[i]));
//   }
//   X[0] *= 3;
// }
// 
// // [[Rcpp::export()]]
// void m3(Map<MatrixXd> X){
//   X *= 3;
// }
// // [[Rcpp::export()]]
// void m3b(Map<MatrixXd> X_){
//   MatrixXd X = X_;
//   Map<MatrixXd> X2(X.data(),X.rows(),X.cols());
//   m3(X2);
// }
// 
// // [[Rcpp::export()]]
// void m4(MatrixXd &X){
//   X *= 3;
//   Rcout << X.coeff(0,0);
// }
// 
// // [[Rcpp::export()]]
// void m4b(Map<MatrixXd> X_){
//   MatrixXd X(X_);
//   m4(X);
//   Rcout << X.coeff(0,0);
// }
// 
// // [[Rcpp::export()]]
// void t(MatrixXd X, VectorXd y){
//   if(y.size() != X.cols()) stop("asdf");
// }

// // [[Rcpp::export()]]
// MatrixXd t1(Map<MatrixXd> R,Map<MatrixXd> Y){
//   return(R.triangularView<Upper>().solve(Y));
// }
// 
// // [[Rcpp::export()]]
// MatrixXd invt(Map<MatrixXd> Rinv,Map<MatrixXd> Y){
//   return(Rinv.triangularView<Upper>() * Y);
// }
// 
// 
// // [[Rcpp::export()]]
// void temp(SEXP Rinv_,Map<MatrixXd> Y, bool sparseR = false){
//   if(Rf_isMatrix(Rinv_) ){
//     Rcout << "yes";
//   } else{
//     Rcout << "no";
//   }
// }
//   
// // [[Rcpp::export()]]
// MatrixXd invt2(SEXP Rinv_,Map<MatrixXd> Y, bool sparseR = false){
//   MatrixXd result;
//   if(sparseR) {
//     MSpMat Rinv = as<MSpMat>(Rinv_);
//     result = Rinv.triangularView<Upper>() * Y;
//   } else{
//     Map<MatrixXd> Rinv = as<Map<MatrixXd> >(Rinv_);
//     result = Rinv.triangularView<Upper>() * Y;
//   }
//   return(result);
// }
// 
// 
// // [[Rcpp::export()]]
// MatrixXd m1(Map<MatrixXd> X, Map<MatrixXd> Y) {
//   MatrixXd result(X.rows(),Y.cols());
//   result = X*Y;
//   return(result);
// }
// // [[Rcpp::export()]]
// MatrixXd m1b(Map<MatrixXd> X, Map<MatrixXd> Y) {
//   MatrixXd result(X.rows(),Y.cols());
//   for(int i = 0; i < Y.cols(); i++){
//     MatrixXd Yi(X.rows(),1);
//     // Yi = Y.col(i);
//     // result.col(i) = X*Y.block(0,i,Y.rows(),1);
//     result.block(0,i,Y.rows(),1) = X*Y.col(i);
//   }
//   return(result);
// }
// // [[Rcpp::export()]]
// MatrixXd m1c(Map<MatrixXd> Xt, Map<MatrixXd> Y) {
//   MatrixXd result(Xt.cols(),Y.cols());
//   for(int i = 0; i < Y.cols(); i++){
//     result.col(i) = Xt.transpose()*Y.col(i);
//   }
//   return(result);
// }
// 
// // [[Rcpp::export()]]
// MatrixXd m2(Map<MatrixXd> X, Map<MatrixXd> Y) {
//   int n = X.cols();
//   int b = Y.rows()/n;
//   int p = Y.cols();
//   MatrixXd result(Y.rows(),Y.cols());
//   for(int i = 0; i < b; i++) {
//     result.block(n*i,0,n,p) = X * Y.block(n*i,0,n,p);
//   }
//   return(result);
// }
// 
// // [[Rcpp::export()]]
// MatrixXd m3(Map<MatrixXd> X, Map<MatrixXd> Y) {
//   int n = X.cols();
//   int b = Y.rows()/n;
//   int p = Y.cols();
//   MatrixXd result(Y.rows(),Y.cols());
//   for(int i = 0; i < p; i++) {
//     MatrixXd Yi = Map<MatrixXd>(Y.col(i).data(),n,b);
//     MatrixXd r = X * Yi;
//     // result.col(i) = Map<VectorXd>(r.data(),n*b);
//   }
//   return(result);
// }
// 
// // [[Rcpp::export()]]
// MatrixXd m4(MatrixXd &X, MatrixXd &Y) {
//   return(X*Y);
// }
// // [[Rcpp::export()]]
// MatrixXd m4b(MatrixXd X, MatrixXd Y) {
//   return(X*Y);
// }
