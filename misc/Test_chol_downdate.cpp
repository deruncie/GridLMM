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

MatrixXd block_cholesky2(MatrixXd L_A, // Cholesky L matrix for A matrix (lower-left block)
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


VectorXd vectorize_chol2(MatrixXd L) {
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




MatrixXd solve_cholesky2(MatrixXd L, MatrixXd X) {
  // solve(LLt,X)
  // L is lower-triangular
  return(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)));
}


struct SS_result{
  MatrixXd beta_hats;
  VectorXd RSSs;
  double V_log_det;
  double V_star_inv_log_det;
  VectorXd V_star_L;
};

SS_result GridLMM_SS(
  MatrixXd &Y_std,
  MatrixXd &X_std,
  MatrixXd &L,
  double V_log_det
) {
  SS_result result;

  result.V_log_det = V_log_det;
  result.V_star_L = vectorize_chol2(L);
  result.V_star_inv_log_det = 2*L.diagonal().array().log().sum();

  MatrixXd b_hat = solve_cholesky2(L,X_std.transpose() * Y_std);
  MatrixXd eta_std = Y_std - X_std * b_hat;

  result.beta_hats = b_hat;
  result.RSSs = eta_std.cwiseProduct(eta_std).colwise().sum();

  return(result);
}


Rcpp::List collect_SS_results(
    std::vector<SS_result> results
) {
  int p = results.size();
  int m = results[0].RSSs.size();
  int b = results[0].beta_hats.size()/m;
  
  MatrixXd RSSs(m,p);
  MatrixXd beta_hats(b*m,p);
  MatrixXd V_star_L(b*(b+1)/2,p);
  VectorXd V_star_inv_log_det(p);
  VectorXd V_log_det(p);
  
  for(int i = 0; i < p; i++){
    SS_result results_i = results[i];
    beta_hats.col(i) = Map<VectorXd>(results_i.beta_hats.data(),results_i.beta_hats.size());
    RSSs.col(i) = results_i.RSSs;
    V_star_L.col(i) = results_i.V_star_L;
    V_star_inv_log_det(i) = results_i.V_star_inv_log_det;
    V_log_det(i) = results_i.V_log_det;
  }
  
  return(Rcpp::List::create(Named("beta_hats") = beta_hats,
                            Named("RSSs") = RSSs,
                            Named("V_log_dets") = V_log_det,
                            Named("V_star_inv_log_det") = V_star_inv_log_det,
                            Named("V_star_L") = V_star_L));
}


// [[Rcpp::export()]]
Rcpp::List GridLMM_SS_dense_c2(
    Map<MatrixXd> Y,
    SEXP chol_Vi_R_,
    Map<MatrixXd> X_cov,
    Map<MatrixXd> X,
    Map<MatrixXd> X_test,
    ArrayXi X_indices,
    VectorXd inv_prior_X,
    double V_log_det
) {
  
  bool sparseR = Rf_isMatrix(chol_Vi_R_);
  
  int b_cov = X_cov.cols();
  int b_x = X_test.cols();
  int b = b_cov + b_x;
  int p = X_indices.size();
  int n = Y.rows();
  std::vector<SS_result> results;

  MatrixXd Y_std, X_cov_std, Xi_std;
  if(sparseR) {
    Map<MatrixXd> chol_Vi_R = as<Map<MatrixXd> >(chol_Vi_R_);
    Y_std = chol_Vi_R.transpose().triangularView<Lower>().solve(Y);
    X_cov_std = chol_Vi_R.transpose().triangularView<Lower>().solve(X_cov);
  } else{
    MSpMat chol_Vi_R = as<MSpMat>(chol_Vi_R_);
    Y_std = chol_Vi_R.transpose().triangularView<Lower>().solve(Y);
    X_cov_std = chol_Vi_R.transpose().triangularView<Lower>().solve(X_cov);
  }
  // MatrixXd Y_std = inv_chol_Vi_transpose.triangularView<Lower>() * Y;
  // MatrixXd X_cov_std = inv_chol_Vi_transpose.triangularView<Lower>() * X_cov;
  
  MatrixXd A = X_cov_std.transpose() * X_cov_std;
  A.diagonal() += inv_prior_X.head(b_cov);
  Eigen::LLT<MatrixXd> llt_of_A(A);
  MatrixXd L_A = llt_of_A.matrixL();


  if(p == 0) {
    // If there are no tests to do, just run the null model
    MatrixXd X_std = X_cov_std;
    SS_result result_1 = GridLMM_SS(Y_std,X_std,L_A,V_log_det);
    results.push_back(result_1);
  } else{
    for(int i = 0; i < p; i++) {
      // for each column of X, form the test matrix by multiplying the marker by the design matrix
      MatrixXd Xi = X.col(X_indices[i]-1).asDiagonal() * X_test;
      
      if(sparseR) {
        Map<MatrixXd> chol_Vi_R = as<Map<MatrixXd> >(chol_Vi_R_);
        Xi_std = chol_Vi_R.transpose().triangularView<Lower>().solve(Xi);
      } else{
        MSpMat chol_Vi_R = as<MSpMat>(chol_Vi_R_);
        Xi_std = chol_Vi_R.transpose().triangularView<Lower>().solve(Xi);
      }
      // MatrixXd Xi_std = inv_chol_Vi_transpose.triangularView<Lower>() * Xi;
      MatrixXd X_std(n,b);
      X_std << X_cov_std,Xi_std;
      
      // calculate L
      MatrixXd B = X_cov_std.transpose() * Xi_std;
      MatrixXd D = Xi_std.transpose() * Xi_std;
      D.diagonal() += inv_prior_X.tail(b_x);
      MatrixXd L = block_cholesky2(L_A, B, D);
      
      // Then solve the model
      SS_result result_i = GridLMM_SS(Y_std,X_std,L,V_log_det);
      results.push_back(result_i);
    }
  }

  return(collect_SS_results(results));
}