#include <RcppEigen.h>


using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::ArrayXd;                  // variable size vector, double precision
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;

using namespace Rcpp;
using namespace Eigen;




// [[Rcpp::export()]]
MatrixXd chol_update2(MatrixXd L, MatrixXd X, int sign) {
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
MatrixXd chol_update2s(MatrixXd L, MatrixXd X, int sign) {
  int n = L.rows();
  Map<MatrixXd> L_new(L.data(),n,n);
  for(int i = 0; i < X.cols(); i++){
    int sum_nzero = 0;
    // Rcout << X << std::endl << std::endl;
    for(int k = 0; k < n; k++){
      double x_k = X.coeffRef(k,i);
      // if(x_k != 0.0) {
      if(std::abs(x_k - 0.0) > 1e-10) {
        // Rcout << i << " " << k << " " << x_k << std::endl;
        sum_nzero++;
        double L_kk = L_new.coeffRef(k,k);
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
    // Rcout << sum_nzero << std::endl;
  }
  return(L_new);
}
// 
// 
// MatrixXd block_cholesky2(MatrixXd L_A, // Cholesky L matrix for A matrix (lower-left block)
//                         MatrixXd B, MatrixXd D){
//   // block cholesky decomposition of a symmetric square matrix: [A B;B^t D]
//   // returns L st L*L^T = original matrix
//   MatrixXd L_A_invB = L_A.triangularView<Lower>().solve(B);
//   MatrixXd Q = D - L_A_invB.transpose() * L_A_invB;
//   Eigen::LLT<MatrixXd> llt_of_Q(Q);
//   MatrixXd QL = llt_of_Q.matrixL();
//   MatrixXd L(L_A.rows() + D.rows(),L_A.rows() + D.rows());
//   MatrixXd z = MatrixXd::Zero(L_A.rows(),D.cols());
//   L << L_A,z,L_A_invB.transpose(),QL;
//   return(L);
// }
// 
// // [[Rcpp::export()]]
// VectorXd log_det_of_XtX2(
//     Map<MatrixXd> X_cov,
//     Rcpp::List X_tests,
//     ArrayXi X_indices   // 1-based index of the tests to actually calculate
// ){
//   
//   int b_x = X_tests.length();
//   int n = X_cov.rows();
//   
//   MatrixXd A = X_cov.transpose() * X_cov;
//   Eigen::LLT<MatrixXd> llt_of_A(A);
//   MatrixXd L_A = llt_of_A.matrixL();
//   
//   std::vector<Map<MatrixXd> > Xs;
//   for(int i = 0; i < b_x; i++){
//     Xs.push_back(as<Map<MatrixXd> >(X_tests[i]));
//     if(Xs[i].cols() != Xs[0].cols()) stop("Different numbers of columns in X_list matrices");
//   }
//   int p = X_indices.size();
//   
//   if(p == 0){
//     VectorXd log_det(1);
//     log_det(0) = 2*L_A.diagonal().array().log().sum();
//     return(log_det);
//   }
//   
//   
//   VectorXd log_det(p);
//   MatrixXd Xi(n,b_x);
//   for(int i = 0; i < p; i++){
//     for(int j = 0; j < b_x; j++){
//       Xi.col(j) = Xs[j].col(X_indices[i]-1);
//     }
//     
//     MatrixXd B = X_cov.transpose() * Xi;
//     MatrixXd D = Xi.transpose() * Xi;
//     MatrixXd L = block_cholesky2(L_A, B, D);
//     
//     log_det(i) = 2*L.diagonal().array().log().sum();
//   }
//   
//   return(log_det);
// }
// 
// void load_list(
//     Rcpp::List X_list_,
//     ArrayXi X_indices
//     ) {
//   
//   int b_x = X_list_.size();
//   // Process X_list
//   std::vector<Map<MatrixXd> > X_list;
//   for(int j = 0; j < b_x; j++) {
//     X_list.push_back(as<Map<MatrixXd> >(X_list_[j]));
//     if(X_list[j].cols() != X_list[0].cols()) stop("Different numbers of columns in X_list matrices");
//     // if(X_list[j].rows() != n) stop("Wrong number of rows in X_list matrices");
//   }
//   int n = X_list[0].rows();
//   int p = X_indices.size();
//   
//   MatrixXd X_std(n,p*b_x);
//   for(int i = 0; i < p; i++) {
//     // reform X_list into a wide matrix of n x b_x design matrices
//     int index_i = X_indices[i]-1;
//     for(int j = 0; j < b_x; j++) {
//       X_std.col(i*b_x + j) = X_list[j].col(index_i);
//     }
//   }
//   
//   // Rcout << X_list[0].col(X_list[0].cols()-1) << std::endl;
//   // Rcout << X_list[0].cols() << std::endl;
// }
// 
// MatrixXd ginv_LDLt(Eigen::LDLT<MatrixXd,Lower> ldlt_K, MatrixXd Y) {
//   if(Y.rows() != ldlt_K.rows()) stop("Wrong dimensions of Y");
//   MatrixXd L = ldlt_K.matrixL();
//   MatrixXd I = MatrixXd::Identity(ldlt_K.rows(), ldlt_K.rows());
//   MatrixXd P = ldlt_K.transpositionsP() * I;
//   VectorXd d = ldlt_K.vectorD();
//   int p = d.size();
//   VectorXd dinv = VectorXd::Zero(p);
//   for(int i = 0; i < p; i++) {
//     if(d[i] > 1e-9) dinv[i] = 1/d[i];
//   }
//   MatrixXd LinvP = L.triangularView<Lower>().solve(P);
//   return(LinvP.transpose() * dinv.asDiagonal() * LinvP*Y);
// }
// 
// // [[Rcpp::export()]]
// MatrixXd solve_LDLt(Map<MatrixXd> K, Map<MatrixXd> y){
//   // MatrixXd K = X.transpose()*X;
//   Eigen::LDLT<MatrixXd,Lower> ldlt_K(K);
//   // MatrixXd Z = ldlt_K.solve(y);
//   MatrixXd Z = ginv_LDLt(ldlt_K,y);
//   return(Z);
// }
// 
// 
// // [[Rcpp::export()]]
// List LDLt_downdate(Map<MatrixXd> K, Map<MatrixXd> X) {
//   Eigen::LDLT<MatrixXd,Lower> ldlt_K(K);
//   int k = X.cols();
//   int p = X.rows();
//   // Rcout << ldlt_K.isPositive() << std::endl;
//   for(int i = 0; i< k; i++) {
//     // int sign = 1;
//     // if(X.col(i).minCoeff() < 0) sign = -1;
//     // ldlt_K.rankUpdate(X.col(i).cwiseAbs().cwiseSqrt(),sign);
//     ldlt_K.rankUpdate(X.col(i));
//     // Rcout << ldlt_K.isPositive() << std::endl;
//   }
//   // Rcout << ldlt_K.isPositive() << std::endl;
//   ldlt_K.setZero();
//   MatrixXd L = ldlt_K.matrixL();
//   VectorXd d = ldlt_K.vectorD();
//   MatrixXd I = MatrixXd::Identity(ldlt_K.rows(), ldlt_K.rows());
//   MatrixXd P = ldlt_K.transpositionsP() * I;
//   MatrixXd K2 = ldlt_K.reconstructedMatrix();
//   return(Rcpp::List::create(Named("L") = L,
//                             Named("d") = d,
//                             Named("P") = P,
//                             Named("K2") = K2
//   ));
// }
// // 
// // // [[Rcpp::export()]]
// // VectorXd Calculate_qt(Map<MatrixXd> X, Map<MatrixXd> beta,Map<VectorXd> lambdas){
// //   // This function calculates tr(X(t(X) %*% X + n*lambda*W^-)^- t(X)), which is the number of effective parameters, part of the calculation of GCV.
// //   // It uses the fact that W^- is diagonal, with most elements zero, so several low-rank updates can be made to LDLt = t(X) %*% X
// //   int n = X.rows();
// //   int p = beta.rows();
// //   if(X.cols() != p) stop("Wrong dimensions of beta or X");
// //   int l = lambdas.size();
// //   if(beta.cols() != l) stop("Wrong length of lambdas");
// //   MatrixXd nlambda_beta_abs =(1.0/n)*(beta*(lambdas.cwiseInverse().asDiagonal())).cwiseAbs();
// //   MatrixXd K = X.transpose() * X;
// //   Eigen::LDLT<MatrixXd,Lower> ldlt_K_base(K);
// //   VectorXd results = VectorXd::Zero(l);
// //   for(int i = 0; i < l; i++) {
// //     Eigen::LDLT<MatrixXd,Lower> ldlt_K = ldlt_K_base;
// //     for(int j = 0; j < p; j++) {
// //       if(nlambda_beta_abs.coeffRef(j,i) > 0) {
// //         VectorXd w = VectorXd::Zero(p);
// //         w[j] = 1.0/std::sqrt(nlambda_beta_abs.coeffRef(j,i));
// //         ldlt_K.rankUpdate(w);
// //       }
// //     }
// //     MatrixXd S = X * ldlt_K.solve(X.transpose());
// //     results[i] = S.diagonal().sum();
// //   }
// //   return(results);
// // }
// 
// 
// // [[Rcpp::export()]]
// List LDLT_base(MatrixXd K) {
//   LDLT<MatrixXd,Lower> ldlt(K);
//   MatrixXd mLDLT = ldlt.matrixLDLT();
//   // Rcout << m << std::endl;
//   // MatrixXd L = ldlt.matrixL();
//   // VectorXd d = ldlt.vectorD();
//   VectorXi t = ldlt.transpositionsP().indices();
//   return(Rcpp::List::create(Named("mLDLT") = mLDLT,
//                     // Named("L") = L,
//                     //  Named("d") = d,
//                      Named("t") = t));
// }
// 
// // [[Rcpp::export()]]
// List LDLT_u1(MatrixXd K,MatrixXd X) {
//   LDLT<MatrixXd,Lower> ldlt(K);
//   int k = X.cols();
//   for(int i = 0; i < k; i++) {
//     ldlt.rankUpdate(X.col(i));
//   }
//   MatrixXd mLDLT = ldlt.matrixLDLT();
//   // Rcout << m << std::endl;
//   // MatrixXd L = ldlt.matrixL();
//   // VectorXd d = ldlt.vectorD();
//   VectorXi t = ldlt.transpositionsP().indices();
//   return(Rcpp::List::create(Named("mLDLT") = mLDLT,
//                             // Named("L") = L,
//                             //  Named("d") = d,
//                             Named("t") = t));
// }
// 
// // [[Rcpp::export()]]
// List LDLT_update(List LDLt_old, MatrixXd X) {
//   int k = X.cols();
//   int n = X.rows();
//   // Rcout << n << std::endl;
//   MatrixXd mLDLT = as<MatrixXd>(LDLt_old["mLDLT"]);
//   MatrixXd L; //= as<MatrixXd>(LDLt_old["L"]);
//   VectorXd d; //= as<VectorXd>(LDLt_old["d"]);
//   VectorXi t = as<VectorXi>(LDLt_old["t"]);
//   // Rcout << L << std::endl;
//   // Rcout << d << std::endl;
//   LDLT<MatrixXd,Lower> ldlt(n);
//   ldlt.matrixLDLT().const_cast_derived() = mLDLT;
//   ldlt.transpositionsP().indices().const_cast_derived() = t;
//   // t = ldlt.transpositionsP().indices();
//   // Rcout << t << std::endl;
//   // d = ldlt.vectorD();
//   // Rcout << d << std::endl;
//   // ldlt.matrixL().nestedExpression().const_cast_derived() = L;
//   // ldlt.vectorD().nestedExpression().const_cast_derived() = d.asDiagonal();
//   // d = ldlt.vectorD();
//   // Rcout << d << std::endl;
//   // L = ldlt.matrixL();
//   // Rcout << L << std::endl;
//   // ldlt.vectorD().fill(d);
//   // ldlt.vectorD().coeff(0) = d(0);
//   // ldlt.vectorD().nestByValue().const_cast_derived() = d;
//   // ldlt.vectorD().nestedExpression().const_cast_derived() = d.asDiagonal();
//   // diagonal().const_cast_derived() = d;
//   // const_cast_derived() = d;
//   // nestedExpression().const_cast_derived() = d.asDiagonal();
//   // ldlt.matrixL().nestedExpression().const_cast_derived() = L;
//   // nestedExpression().const_cast_derived() = d;
//   // L = ldlt.matrixL();
//   // Rcout << L << std::endl;
//   // d = ldlt.vectorD();
//   // Rcout << d << std::endl;
//   // t = ldlt.transpositionsP().indices();
//   // Rcout << t << std::endl;
//   // L = ldlt.matrixL();
//   // Rcout << L << std::endl;
//   
//   // for(int i = 0; i < k; i++) {
//   //   ldlt.rankUpdate(X.col(i));
//   // }
//   L = ldlt.matrixL();
//   d = ldlt.vectorD();
//   t = ldlt.transpositionsP().indices();
//   return(Rcpp::List::create(Named("L") = L,
//                             Named("d") = d,
//                             Named("t") = t));
// }
// 
// 
// // [[Rcpp::export()]]
// void asdf(VectorXi indices) {
//   PermutationMatrix<Dynamic> P(indices);
//   MatrixXd I = MatrixXd::Identity(indices.size(),indices.size());
//   MatrixXd M = P*I;
//   Rcout << M << std::endl;
//   Rcout << P.indices() << std::endl;
// }
// 
// // [[Rcpp::export()]]
// MatrixXd LDLt2(MatrixXd X) {
//   LDLT<MatrixXd,Lower> ldlt(X);
//   // ldlt.compute();
//   MatrixXd I = MatrixXd::Identity(ldlt.rows(), ldlt.rows());
//   MatrixXd P = ldlt.transpositionsP() * I;
//   Rcout << P << std::endl;
//   Rcout << ldlt.transpositionsP().coeff(0) << std::endl;
//   Rcout << ldlt.transpositionsP().coeff(1) << std::endl;
//   Rcout << ldlt.transpositionsP().coeff(2) << std::endl;
//   Eigen::Transpositions<Dynamic> P2(3); 
//   VectorXi p(3);
//   p <<  2,1,2;
//   P2.indices() = p;
//   ldlt.transpositionsP().indices().const_cast_derived() = p;
//   Rcout << ldlt.transpositionsP().coeff(0) << std::endl;
//   Rcout << ldlt.transpositionsP().coeff(1) << std::endl;
//   Rcout << ldlt.transpositionsP().coeff(2) << std::endl;
//     // const_cast_derived() = P2;
//   P = ldlt.transpositionsP() * I;
//   Rcout << P << std::endl;
//   // ldlt.transpositionsP().nestedExpression().
//   // ldlt.transpositionsP().coeff(1) = 1;
//   // ldlt.transpositionsP().coeff(2) = 1;
//   // ldlt.matrixU().nestedExpression().const_cast_derived()
//   MatrixXd L = ldlt.matrixL();
//   VectorXd d = ldlt.vectorD();
//   // VectorXd p = ldlt.transpositionsP().indices().cast(double);
//   // Rcout << p << std::endl;
//   // Rcout << ldlt.vectorD() << std::endl;
//   // Rcout << ldlt.transpositionsP() << std::endl;
//   return(ldlt.matrixL());
// }
// 
// // // [[Rcpp::export()]]
// // List LKLt_update(List LDLt_old,MatrixXd X) {
// //   MatrixXd L = as<MatrixXd>(LDLt_old["L"]);
// //   MSpMat P = as<MSpMat>(LDLt_old["P"]);
// //   VectorXd d = as<VectorXd>(LDLt_old["d"]);
// //   int n = X.rows();
// //   LDLT<MatrixXd,Lower> ldlt(n);
// //   ldlt.matrixL().nestedExpression().const_cast_derived() = L;
// //   ldlt.vectorD().nestedExpression().const_cast_derived() = d;
// //   // ldlt.transpositionsP().indices().const_cast_derived()
// // }
// 
// // [[Rcpp::export()]]
// MatrixXd my_cholUpdate(Map<MatrixXd> R, MatrixXd X) {
//   int n = X.rows();
//   LLT<MatrixXd,Lower> llt(n);
//   // llt.matrixL().nestedExpression().const_cast_derived() = VectorXd::Ones(n).asDiagonal();
//   llt.matrixL().nestedExpression().const_cast_derived() = R.transpose();
//   for(int i = 0; i < X.cols(); i++) {
//     llt.rankUpdate(X.col(i));
//   }
//   MatrixXd L = llt.matrixL();
//   return(L);
// }
// 
// // // [[Rcpp::export()]]
// // MatrixXd usemyChol1(Map<MatrixXd> R) {
// //   int n = R.rows();
// //   LLT<Ref<MatrixXd>,Lower> llt(n);
// //   llt.matrixL().nestedExpression().const_cast_derived() = R.transpose();
// //   MatrixXd L = llt.matrixL();
// //   return(L);
// // }
// // // [[Rcpp::export()]]
// // MatrixXd usemyChol2(Map<MatrixXd> R) {
// //   // int n = R.rows();
// //   LLT<MatrixXd,Lower> llt(R);
// //   llt.matrixL().nestedExpression().const_cast_derived() = R.transpose();
// //   MatrixXd L = llt.matrixL();
// //   return(L);
// // }
// // // [[Rcpp::export()]]
// // void usemyChol3(MatrixXd R) {
// //   // int n = R.rows();
// //   // LLT<MatrixXd,Lower> llt(n);
// //   // llt.matrixL().nestedExpression().const_cast_derived() = R.transpose();
// //   // MatrixXd L = llt.matrixL();
// //   // MatrixXd L = R;
// //   // return(L);
// // }
// // // [[Rcpp::export()]]
// // void usemyChol4(Map<MatrixXd> R) {
// //   // R.triangularView<Lower>().nestedExpression().
// //   // int n = R.rows();
// //   // LLT<MatrixXd,Lower> llt(n);
// //   // llt.matrixL().nestedExpression().const_cast_derived() = R.transpose();
// //   // MatrixXd L = llt.matrixL();
// //   // Map<MatrixXd> L = R;
// //   // return(R);
// // }
// 
// 
// // void sdi(Rcpp::List A, IntegerVector b){
// //   IntegerVector a(0);
// //   for(int i = 0; i < A.size(); i++) {
// //     Rcout << Rf_isInteger(A[i]) << std::endl;
// //     if(Rf_isInteger(A[i])) {
// //       a = as<IntegerVector>(A[i]);
// //     }
// //   }
// //   Rcout << a << std::endl;
// //   // Rcout << Rf_isInteger(a_) << std::endl;
// //   // if(Rf_isInteger(a_)){
// //   //   IntegerVector a = as<IntegerVector>(a_);
// //   //   return(setdiff(a,b));
// //   // }
// //   // return(b);
// // }
// 
// 
// // void chol_update_R_inplace2(MatrixXd &R, MatrixXd X, VectorXd weights) {
// //   int n = R.rows();
// //   if(X.rows() != n) stop("Wrong dimension of X for downdating R");
// //   if(weights.size() != X.cols()) stop("wrong length of weights for downdating R");
// //   for(int i = 0; i < X.cols(); i++){
// //     VectorXd Xi = X.col(i);
// //     double weights_i = weights[i];
// //     for(int k = 0; k < n; k++){
// //       double R_kk = R.coeffRef(k,k);
// //       double x_k = Xi.coeffRef(k);
// //       double r = sqrt(R_kk*R_kk + weights_i * x_k*x_k);
// //       double c = r / R_kk;
// //       double s = x_k / R_kk;
// //       R.coeffRef(k,k) = r;
// //       if(k < (n-1)) {
// //         for(int j=k+1; j<n; j++) {
// //           double R_kj = R.coeffRef(k,j);
// //           double X_ij = Xi.coeff(j);
// //           R_kj += weights_i * s * X_ij;
// //           R_kj /= c;
// //           X_ij *= c;
// //           X_ij -= s*R_kj;
// //           R.coeffRef(k,j) = R_kj;
// //           Xi.coeffRef(j) = X_ij;
// //         }
// //         // R.block(k,k+1,1,n-k-1) = (R.block(k,k+1,1,n-k-1) + weights_i * s * Xi.tail(n-k-1).transpose())/c;
// //         // Xi.tail(n-k-1) = c*Xi.tail(n-k-1) - s*R.block(k,k+1,1,n-k-1).transpose();
// //         // R.block(k,k+1,1,n-k-1) += weights_i * s * Xi.tail(n-k-1).transpose();
// //         // R.block(k,k+1,1,n-k-1) /= c;
// //         // Xi.tail(n-k-1) *= c;
// //         // Xi.tail(n-k-1) -= s*R.block(k,k+1,1,n-k-1).transpose();
// //         }
// //     }
// //   }
// // }
// // 
// // 
// // // [[Rcpp::export()]]
// // MatrixXd chol_update_L2(MatrixXd L, MatrixXd X, VectorXd weights) {
// //   MatrixXd R = L.transpose();
// //   chol_update_R_inplace2(R,X,weights);
// //   return(R.transpose());
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd chol_update_R2(MatrixXd R, MatrixXd X, VectorXd weights) {
// //   chol_update_R_inplace2(R,X,weights);
// //   return(R);
// // }
// // 
// // void chol_update_L_inplace2(MatrixXd &L, MatrixXd X, VectorXd weights) {
// //   int n = L.rows();
// //   if(X.rows() != n) stop("Wrong dimension of X for downdating L");
// //   if(weights.size() != X.cols()) stop("wrong length of weights for downdating L");
// //   for(int i = 0; i < X.cols(); i++){
// //     VectorXd Xi = X.col(i);
// //     double weights_i = weights[i];
// //     for(int k = 0; k < n; k++){
// //       double L_kk = L.coeffRef(k,k);
// //       double x_k = Xi.coeffRef(k);
// //       double r = sqrt(L_kk*L_kk + weights_i * x_k*x_k);
// //       double c = r / L_kk;
// //       double s = x_k / L_kk;
// //       L.coeffRef(k,k) = r;
// //       if(k < (n-1)) {
// //         for(int j=k+1; j<n; j++) {
// //           double L_jk = L.coeffRef(j,k);
// //           double X_ij = Xi.coeff(j);
// //           L_jk += weights_i * s * X_ij;
// //           L_jk /= c;
// //           X_ij *= c;
// //           X_ij -= s*L_jk;
// //           L.coeffRef(j,k) = L_jk;
// //           Xi.coeffRef(j) = X_ij;
// //         }
// //         // L.block(k+1,k,n-k-1,1) = (L.block(k+1,k,n-k-1,1) + weights_i * s * Xi.tail(n-k-1))/c;
// //         // Xi.tail(n-k-1) = c*Xi.tail(n-k-1) - s*L.block(k+1,k,n-k-1,1);
// //         // R.block(k,k+1,1,n-k-1) += weights_i * s * Xi.tail(n-k-1).transpose();
// //         // R.block(k,k+1,1,n-k-1) /= c;
// //         // Xi.tail(n-k-1) *= c;
// //         // Xi.tail(n-k-1) -= s*R.block(k,k+1,1,n-k-1).transpose();
// //       }
// //     }
// //   }
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd chol_update_L3(MatrixXd L, MatrixXd X, VectorXd weights) {
// //   chol_update_L_inplace2(L,X,weights);
// //   return(L);
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd chol_solve1(Map<MatrixXd> chol_R, Map<MatrixXd> Y){
// //   return(chol_R.transpose().triangularView<Lower>().solve(Y));
// // }
// // // [[Rcpp::export()]]
// // MatrixXd chol_solve2(Map<MatrixXd> chol_L, Map<MatrixXd> Y){
// //   return(chol_L.triangularView<Lower>().solve(Y));
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd chol_u1(Map<MatrixXd> A, VectorXd X, int times){
// //   LLT<MatrixXd> lA(A);
// //   MatrixXd L = lA.matrixL();
// //   // VectorXd weights = VectorXd::Constant(1,1.0);
// //   // for(int i = 0; i < times; i++){
// //   //   chol_update_L_inplace2(L, X, weights);
// //   // }
// //   MatrixXd XX = X.replicate(1,times);
// //   VectorXd weights = VectorXd::Constant(times,1.0);
// //   chol_update_L_inplace2(L,XX,weights);
// //   return(L);
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd chol_u2(Map<MatrixXd> A, VectorXd X, int times){
// //   LLT<MatrixXd> lA(A);
// //   VectorXd weights = VectorXd::Constant(1,1.0);
// //   for(int i = 0; i < times; i++){
// //     lA.rankUpdate(X,1.0);
// //   }
// //   MatrixXd L = lA.matrixL();
// //   return(L);
// // }
// // 
// // Rcpp::List chol_cT(Map<MatrixXd> A){
// //   LLT<MatrixXd> lA(A);
// //   return(Rcpp::List::create(lA));
// // }
// // // [[Rcpp::export()]]
// // MatrixXd chol_cL(Rcpp::List L){
// //   LLT<MatrixXd> lA = as<LLT<MatrixXd>>(L[0]);
// //   return(lA.matrixL());
// // }
// 
// 
// // void a1(MatrixXd &X){
// //   X.coeffRef(0,0) += 1;
// //   Rcout << X.coeff(0,0) << std::endl;
// // }
// // 
// // // [[Rcpp::export()]]
// // void a0(MatrixXd X){
// //   Rcout << X.coeff(0,0) << std::endl;
// //   a1(X);
// //   Rcout << X.coeff(0,0) << std::endl;
// //   MatrixXd Y(X);
// //   a1(Y);
// //   Rcout << X.coeff(0,0) << std::endl;
// // }
// // 
// // struct x{
// //   double A=0;
// // };
// // 
// // x a3a(double A){
// //   x X;
// //   X.A = A;
// //   return(X);
// // }
// // void a3(x &X){
// //   // X.A = 3;
// //   X = a3a(5);
// // }
// // 
// // // [[Rcpp::export()]]
// // void a4(){
// //   x X;
// //   Rcout << X.A << std::endl;
// //   a3(X);
// //   Rcout << X.A << std::endl;
// // }
// 
// // // [[Rcpp::export()]]
// // Rcpp::List svd_c(Map<MatrixXd> X) {
// //   Eigen::BDCSVD<MatrixXd> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
// //   // bcdsolve.computeU();
// //   // MatrixXd U = bcdsolve.matrixU();
// //   return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
// //                      Named("d") = bcdsolve.singularValues()));
// // }
// // 
// // 
// // // [[Rcpp::export()]]
// // Rcpp::List svd_c3(SEXP X) {
// //   if(Rf_isMatrix(X)) {
// //     Eigen::BDCSVD<MatrixXd> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
// //     // bcdsolve.computeU();
// //     // MatrixXd U = bcdsolve.matrixU();
// //     return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
// //                               Named("d") = bcdsolve.singularValues()));
// //   } else{
// //     Eigen::BDCSVD<MSpMat> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
// //     // bcdsolve.computeU();
// //     // MatrixXd U = bcdsolve.matrixU();
// //     return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
// //                               Named("d") = bcdsolve.singularValues()));
// //   }
// // }
// // 
// // // [[Rcpp::export()]]
// // Rcpp::List svd_c2(MatrixXd X) {
// //   Eigen::BDCSVD<MatrixXd> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
// //   // bcdsolve.computeU();
// //   // MatrixXd U = bcdsolve.matrixU();
// //   return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
// //                             Named("d") = bcdsolve.singularValues()));
// // }
// // 
// // 
// // // [[Rcpp::export()]]
// // void m1(Map<MatrixXd> X){
// //   X *= 3;
// // }
// // 
// // // [[Rcpp::export()]]
// // void m1b(Map<MatrixXd> X){
// //   MatrixXd X2(X);
// //   X2 *= 3;
// // }
// // 
// // 
// // // [[Rcpp::export()]]
// // void test_list(Rcpp::List X_){
// //   std::vector<Map<MatrixXd> > X;
// //   for(int i = 0; i < X_.length(); i++){
// //     X.push_back(as<Map<MatrixXd> >(X_[i]));
// //   }
// //   X[0] *= 3;
// // }
// // // [[Rcpp::export()]]
// // void test_list2(Rcpp::List X_){
// //   std::vector<MatrixXd> X;
// //   for(int i = 0; i < X_.length(); i++){
// //     X.push_back(MatrixXd());
// //     MatrixXd Xi = as<MatrixXd>(X_[i]);
// //     X.back().swap(as<Map<MatrixXd> >(X_[i]));
// //   }
// //   X[0] *= 3;
// // }
// // 
// // // [[Rcpp::export()]]
// // void m3(Map<MatrixXd> X){
// //   X *= 3;
// // }
// // // [[Rcpp::export()]]
// // void m3b(Map<MatrixXd> X_){
// //   MatrixXd X = X_;
// //   Map<MatrixXd> X2(X.data(),X.rows(),X.cols());
// //   m3(X2);
// // }
// // 
// // // [[Rcpp::export()]]
// // void m4(MatrixXd &X){
// //   X *= 3;
// //   Rcout << X.coeff(0,0);
// // }
// // 
// // // [[Rcpp::export()]]
// // void m4b(Map<MatrixXd> X_){
// //   MatrixXd X(X_);
// //   m4(X);
// //   Rcout << X.coeff(0,0);
// // }
// // 
// // // [[Rcpp::export()]]
// // void t(MatrixXd X, VectorXd y){
// //   if(y.size() != X.cols()) stop("asdf");
// // }
// 
// // // [[Rcpp::export()]]
// // MatrixXd t1(Map<MatrixXd> R,Map<MatrixXd> Y){
// //   return(R.triangularView<Upper>().solve(Y));
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd invt(Map<MatrixXd> Rinv,Map<MatrixXd> Y){
// //   return(Rinv.triangularView<Upper>() * Y);
// // }
// // 
// // 
// // // [[Rcpp::export()]]
// // void temp(SEXP Rinv_,Map<MatrixXd> Y, bool sparseR = false){
// //   if(Rf_isMatrix(Rinv_) ){
// //     Rcout << "yes";
// //   } else{
// //     Rcout << "no";
// //   }
// // }
// //   
// // // [[Rcpp::export()]]
// // MatrixXd invt2(SEXP Rinv_,Map<MatrixXd> Y, bool sparseR = false){
// //   MatrixXd result;
// //   if(sparseR) {
// //     MSpMat Rinv = as<MSpMat>(Rinv_);
// //     result = Rinv.triangularView<Upper>() * Y;
// //   } else{
// //     Map<MatrixXd> Rinv = as<Map<MatrixXd> >(Rinv_);
// //     result = Rinv.triangularView<Upper>() * Y;
// //   }
// //   return(result);
// // }
// // 
// // 
// // // [[Rcpp::export()]]
// // MatrixXd m1(Map<MatrixXd> X, Map<MatrixXd> Y) {
// //   MatrixXd result(X.rows(),Y.cols());
// //   result = X*Y;
// //   return(result);
// // }
// // // [[Rcpp::export()]]
// // MatrixXd m1b(Map<MatrixXd> X, Map<MatrixXd> Y) {
// //   MatrixXd result(X.rows(),Y.cols());
// //   for(int i = 0; i < Y.cols(); i++){
// //     MatrixXd Yi(X.rows(),1);
// //     // Yi = Y.col(i);
// //     // result.col(i) = X*Y.block(0,i,Y.rows(),1);
// //     result.block(0,i,Y.rows(),1) = X*Y.col(i);
// //   }
// //   return(result);
// // }
// // // [[Rcpp::export()]]
// // MatrixXd m1c(Map<MatrixXd> Xt, Map<MatrixXd> Y) {
// //   MatrixXd result(Xt.cols(),Y.cols());
// //   for(int i = 0; i < Y.cols(); i++){
// //     result.col(i) = Xt.transpose()*Y.col(i);
// //   }
// //   return(result);
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd m2(Map<MatrixXd> X, Map<MatrixXd> Y) {
// //   int n = X.cols();
// //   int b = Y.rows()/n;
// //   int p = Y.cols();
// //   MatrixXd result(Y.rows(),Y.cols());
// //   for(int i = 0; i < b; i++) {
// //     result.block(n*i,0,n,p) = X * Y.block(n*i,0,n,p);
// //   }
// //   return(result);
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd m3(Map<MatrixXd> X, Map<MatrixXd> Y) {
// //   int n = X.cols();
// //   int b = Y.rows()/n;
// //   int p = Y.cols();
// //   MatrixXd result(Y.rows(),Y.cols());
// //   for(int i = 0; i < p; i++) {
// //     MatrixXd Yi = Map<MatrixXd>(Y.col(i).data(),n,b);
// //     MatrixXd r = X * Yi;
// //     // result.col(i) = Map<VectorXd>(r.data(),n*b);
// //   }
// //   return(result);
// // }
// // 
// // // [[Rcpp::export()]]
// // MatrixXd m4(MatrixXd &X, MatrixXd &Y) {
// //   return(X*Y);
// // }
// // // [[Rcpp::export()]]
// // MatrixXd m4b(MatrixXd X, MatrixXd Y) {
// //   return(X*Y);
// // }
