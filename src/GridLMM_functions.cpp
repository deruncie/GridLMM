#include <math.h>
#include "GridLMM_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// [[Rcpp::export()]]
Rcpp::List premultiply_list_of_matrices(MSpMat Qt, Rcpp::List X_list){
  int p = X_list.length();
  Rcpp::List X_list_new;
  for(int i = 0; i < p; i++){
    MatrixXd Xi = Rcpp::as<Map<MatrixXd> >(X_list[i]);
    X_list_new.push_back(Qt * Xi);
  }
  return(X_list_new);
}  

double log_det_from_LDLt(Eigen::LDLT<MatrixXd> ldlt){
  double log_det = 0;
  ArrayXd D = ldlt.vectorD();
  ArrayXd D_abs = D.abs();
  for(int i = 0; i < D_abs.size(); i++){
    if(D_abs[i] > 0) log_det += log(D_abs[i]);
  }
  return(log_det);
}


MatrixXd block_cholesky(MatrixXd L_A, // Cholesky L matrix for A matrix (upper-left block)
                        MatrixXd B, MatrixXd D){
  // block cholesky decomposition of a symmetric square matrix: [A B;B^t D]
  // returns L st L*L^T = original matrix
  MatrixXd L_A_invB = L_A.triangularView<Lower>().solve(B);
  MatrixXd Q = D - L_A_invB.transpose() * L_A_invB;
  Eigen::LLT<MatrixXd> llt_of_Q(Q);
  MatrixXd QL = llt_of_Q.matrixL();
  MatrixXd L(L_A.rows() + D.rows(),L_A.rows() + D.rows());
  MatrixXd z = MatrixXd::Zero(L_A.rows(),D.cols());
  L << L_A,z,L_A_invB.transpose(),QL;
  return(L);
}

MatrixXd solve_cholesky(MatrixXd L, MatrixXd X) {
  // solve(LLt,X)
  // L is lower-triangular
  return(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)));
}


VectorXd vectorize_chol(MatrixXd L) {
  int p = L.cols();
  VectorXd vec_chol(p*(p+1)/2);
  int counter = 0;
  for(int i = 0; i < p; i++){
    for(int j = i; j < p; j++){
      vec_chol(counter) = L.coeffRef(j,i);
      counter += 1;
    }
  }
  return(vec_chol);
}

MatrixXd fill_chol_L(VectorXd vec_chol,int p) {
  MatrixXd L = MatrixXd::Zero(p,p);
  int counter = 0;
  for(int i = 0; i < p; i++){
    for(int j = i; j < p; j++){
      L.coeffRef(j,i) = vec_chol(counter);
      counter += 1;
    }
  }
  return(L);
}


// [[Rcpp::export()]]
MatrixXd F_hats(
    Map<MatrixXd> beta_hats,
    Map<MatrixXd> RSSs,
    Map<MatrixXd> V_star_L,
    int n,
    int b,
    int m
) {
  int p = beta_hats.cols();
  MatrixXd F_hats(b*m,p);
  
  MatrixXd I = VectorXd::Ones(b).asDiagonal();
  for(int i = 0; i < p; i++){
    MatrixXd L = fill_chol_L(V_star_L.col(i),b);
    ArrayXXd b_cov = solve_cholesky(L,I).diagonal().array().replicate(1,m) * RSSs.col(i).transpose().array().replicate(b,1)/(n-b);
    Map<MatrixXd> b_hat(beta_hats.col(i).data(),b,m);
    MatrixXd F_hat = b_hat.array().pow(2) / b_cov;
    F_hats.col(i) = Map<VectorXd>(F_hat.data(),F_hat.size());
  }
  return(F_hats);
}

// [[Rcpp::export()]]
VectorXd log_det_of_XtX(
    Map<MatrixXd> X_cov,
    Rcpp::List X_tests,
    ArrayXi X_indices   // 1-based index of the tests to actually calculate
){
  
  int b_x = X_tests.length();
  int n = X_cov.rows();
  
  MatrixXd A = X_cov.transpose() * X_cov;
  Eigen::LLT<MatrixXd> llt_of_A(A);
  MatrixXd L_A = llt_of_A.matrixL();
  
  std::vector<MatrixXd> Xs;
  for(int i = 0; i < b_x; i++){
    MatrixXd Xi = Rcpp::as<Map<MatrixXd> >(X_tests[i]);
    Xs.push_back(Xi);
  }
  int p = X_indices.size();
  
  if(p == 0){
    VectorXd log_det(1);
    log_det(0) = 2*L_A.diagonal().array().log().sum();
    return(log_det);
  }
  
  
  VectorXd log_det(p);
  MatrixXd Xi(n,b_x);
  for(int i = 0; i < p; i++){
    for(int j = 0; j < b_x; j++){
      Xi.col(j) = Xs[j].col(X_indices[i]-1);
    }
    
    MatrixXd B = X_cov.transpose() * Xi;
    MatrixXd D = Xi.transpose() * Xi;
    MatrixXd L = block_cholesky(L_A, B, D);
    
    log_det(i) = 2*L.diagonal().array().log().sum();
  }
  
  return(log_det);
}



Rcpp::List GridLMM_SS(
    MatrixXd &Y_std,
    MatrixXd &X_cov_std,
    std::vector<MatrixXd> &X_stds,
    VectorXd &inv_prior_X,  // diagonal prior precision on X_cov and Xs - same for all tests,
    double V_log_det
) {
  
  int b_cov = X_cov_std.cols();
  int b_x = X_stds.size();
  int b = b_cov+b_x;
  int n = Y_std.rows();
  int m = Y_std.cols();
  int p = X_stds[0].cols();
  
  MatrixXd A = X_cov_std.transpose() * X_cov_std;
  A.diagonal() += inv_prior_X.head(b_cov);
  Eigen::LLT<MatrixXd> llt_of_A(A);
  MatrixXd L_A = llt_of_A.matrixL();
  
  if(p == 0) {
    VectorXd V_star_L = vectorize_chol(L_A);
    VectorXd V_star_inv_log_det(1);
    V_star_inv_log_det(0) = 2*L_A.diagonal().array().log().sum();
    MatrixXd b_hat = solve_cholesky(L_A,X_cov_std.transpose() * Y_std);
    MatrixXd beta_hats(Map<VectorXd>(b_hat.data(),b_hat.size(),1));
    
    MatrixXd eta_std = Y_std - X_cov_std * b_hat;
    
    VectorXd RSSs = eta_std.cwiseProduct(eta_std).colwise().sum();
    
    return(Rcpp::List::create(Named("beta_hats") = beta_hats,
                              Named("RSSs") = RSSs,
                              Named("V_log_dets") = V_log_det,
                              Named("V_star_inv_log_det") = V_star_inv_log_det,
                              Named("V_star_L") = V_star_L));
  }
  
  MatrixXd RSSs(m,p);
  MatrixXd beta_hats((b_cov+b_x)*m,p);
  MatrixXd V_star_L(b*(b+1)/2,p);
  VectorXd V_star_inv_log_det(p);
  
  MatrixXd Xi_std(n,b_x);
  for(int i = 0; i < p; i++){
    for(int j = 0; j < b_x; j++){
      Xi_std.col(j) = X_stds[j].col(i);
    }
    
    MatrixXd B = X_cov_std.transpose() * Xi_std;
    MatrixXd D = Xi_std.transpose() * Xi_std;
    D.diagonal() += inv_prior_X.tail(b_x);
    MatrixXd L = block_cholesky(L_A, B, D);
    
    V_star_L.col(i) = vectorize_chol(L);
    V_star_inv_log_det(i) = 2*L.diagonal().array().log().sum();
    
    MatrixXd X_std(n,b);
    X_std << X_cov_std,Xi_std;
    MatrixXd b_hat = solve_cholesky(L,X_std.transpose() * Y_std);
    beta_hats.col(i) = Map<VectorXd>(b_hat.data(),b_hat.size());
    
    MatrixXd eta_std = Y_std - X_std * b_hat;
    
    RSSs.col(i) = eta_std.cwiseProduct(eta_std).colwise().sum();
  }
  
  return(Rcpp::List::create(Named("beta_hats") = beta_hats,
                            Named("RSSs") = RSSs,
                            Named("V_log_dets") = VectorXd::Constant(p,V_log_det),
                            Named("V_star_inv_log_det") = V_star_inv_log_det,
                            Named("V_star_L") = V_star_L));
}

// [[Rcpp::export()]]
Rcpp::List GridLMM_SS_dense_c(
    Map<MatrixXd> Y,
    Map<MatrixXd> inv_chol_Vi_transpose,
    Map<MatrixXd> X_cov,
    Rcpp::List X_tests,
    ArrayXi X_indices,
    VectorXd inv_prior_X,
    double V_log_det
) {
  
  MatrixXd Y_std = inv_chol_Vi_transpose.triangularView<Lower>() * Y;
  
  int b_x = X_tests.length();
  int p = X_indices.size();
  int n = Y.rows();
  std::vector<MatrixXd> X_stds;
  for(int i = 0; i < b_x; i++){
    MatrixXd Xi_full = Rcpp::as<Map<MatrixXd> >(X_tests[i]);
    MatrixXd Xi(n,p);
    for(int j = 0; j < p; j++) {
      Xi.col(j) = Xi_full.col(X_indices[j]-1);
    }
    X_stds.push_back(inv_chol_Vi_transpose.triangularView<Lower>() * Xi);
  }
  
  MatrixXd X_cov_std = inv_chol_Vi_transpose.triangularView<Lower>() * X_cov;
  
  return(GridLMM_SS(Y_std,X_cov_std,X_stds,inv_prior_X,V_log_det));
}

// [[Rcpp::export()]]
Rcpp::List GridLMM_SS_sparse_c(
    Map<MatrixXd> Y,
    MSpMat inv_chol_Vi_transpose,
    Map<MatrixXd> X_cov,
    Rcpp::List X_tests,
    ArrayXi X_indices,
    VectorXd inv_prior_X,
    double V_log_det
) {
  
  MatrixXd Y_std = inv_chol_Vi_transpose.triangularView<Lower>() * Y;
  
  int b_x = X_tests.length();
  int p = X_indices.size();
  int n = Y.rows();
  std::vector<MatrixXd> X_stds;
  for(int i = 0; i < b_x; i++){
    MatrixXd Xi_full = Rcpp::as<Map<MatrixXd> >(X_tests[i]);
    MatrixXd Xi(n,p);
    for(int j = 0; j < p; j++) {
      Xi.col(j) = Xi_full.col(X_indices[j]-1);
    }
    X_stds.push_back(inv_chol_Vi_transpose.triangularView<Lower>() * Xi);
  }
  
  MatrixXd X_cov_std = inv_chol_Vi_transpose.triangularView<Lower>() * X_cov;
  
  return(GridLMM_SS(Y_std,X_cov_std,X_stds,inv_prior_X,V_log_det));
}



// [[Rcpp::export()]]
Rcpp::List GridLMM_SS_downdate(
    Map<MatrixXd> Y,
    Map<MatrixXd> V_inv,
    Map<MatrixXd> X_cov,
    Rcpp::List X_tests,
    VectorXd inv_prior_X, // diagonal prior precision on X_cov and Xs - same for all tests
    Rcpp::List downdate_weights,
    Rcpp::List downdate_Xs,
    ArrayXi X_indices,
    double V_log_det
) {
  
  
  int p = X_indices.size();
  int n = Y.rows();
  int m = Y.cols();
  int b_cov = X_cov.cols();
  // X_tests is a list of n x p matrices. Tests are constructred by pasting together corresponding columns
  // If this is the null model, X_tests will be empty, b_x will be zero, but downdate_Xs will still have an entry per test
  int b_x = X_tests.length();  
  std::vector<MatrixXd> Xs;
  for(int i = 0; i < b_x; i++){
    MatrixXd Xi = Rcpp::as<Map<MatrixXd> >(X_tests[i]);
    Xs.push_back(Xi);
  }
  if(b_x == 1 && Xs[0].cols() == 0) b_x = 0;  // ensure b_x == 0 for null model
  
  int b = b_cov+b_x;
  
  MatrixXd RSSs(m,p);
  MatrixXd beta_hats((b_cov+b_x)*m,p);
  MatrixXd V_star_L(b*(b+1)/2,p);
  VectorXd V_star_inv_log_det(p);
  VectorXd V_log_dets(p);
  
  for(int i = 0; i < p; i++){
    // construct Xi for this test
    MatrixXd Xi(n,b_x);
    for(int j = 0;j < b_x; j++){
      Xi.col(j) = Xs[j].col(X_indices[i]-1);
    }
    
    MatrixXd V_inv_i = V_inv;
    V_log_dets(i) = V_log_det;
    
    VectorXd downdate_weights_i = as<VectorXd>(downdate_weights[X_indices[i]-1]);
    
    if(downdate_weights_i.size() > 0) {
      // do downdate of V_inv based on downdate_Xi
      MatrixXd dXi = Rcpp::as<Map<MatrixXd> >(downdate_Xs[X_indices[i]-1]);
      MatrixXd V_inv_dXi = V_inv * dXi;
      MatrixXd V_inv_i_inner = V_inv_dXi.transpose() * dXi;
      V_inv_i_inner.diagonal() -= downdate_weights_i.cwiseInverse();
      Eigen::LDLT<MatrixXd> ldlt_of_V_inv_i_inner(V_inv_i_inner);
      V_inv_i -= V_inv_dXi * ldlt_of_V_inv_i_inner.solve(V_inv_dXi.transpose());
      V_log_dets(i) += log_det_from_LDLt(ldlt_of_V_inv_i_inner) + downdate_weights_i.array().log().sum();
    }
    
    MatrixXd Y_std = V_inv_i * Y;
    MatrixXd X(n,X_cov.cols() + Xi.cols());
    X << X_cov,Xi;
    MatrixXd X_std = V_inv_i * X;
    
    MatrixXd V_star = X.transpose() * X_std;
    V_star.diagonal() += inv_prior_X;  
    Eigen::LLT<MatrixXd> llt_of_V_star(V_star);
    MatrixXd L = llt_of_V_star.matrixL();
    V_star_L.col(i) = vectorize_chol(L);
    V_star_inv_log_det(i) = 2*L.diagonal().array().log().sum();
    MatrixXd b_hat = solve_cholesky(L,X.transpose() * Y_std);
    beta_hats.col(i) = Map<VectorXd>(b_hat.data(),b_hat.size());
    
    MatrixXd eta = Y - X * b_hat;
    MatrixXd eta_std = Y_std - X_std * b_hat;
    
    for(int j = 0; j < m; j++){
      RSSs.coeffRef(j,i) = (eta.col(j).transpose() * eta_std.col(j)).coeff(0);
    }
  }
  
  return(Rcpp::List::create(Named("beta_hats") = beta_hats,
                            Named("RSSs") = RSSs,
                            Named("V_log_dets") = V_log_dets,
                            Named("V_star_inv_log_det") = V_star_inv_log_det,
                            Named("V_star_L") = V_star_L));
}
