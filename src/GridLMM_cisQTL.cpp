#include "GridLMM_types.h"
using namespace Rcpp;



// [[Rcpp::export()]]
Rcpp::List GridLMM_SS_cisQTL(
    Map<MatrixXd> Y,
    SEXP chol_Vi_R_,
    MatrixXd X_cov,
    Map<MatrixXd> X,
    MatrixXd X_indices,  // 3 columns: 1) column index of Y, 2) starting column index of X, 3) # of columns of X to test
    Rcpp::List downdate_Xs = Rcpp::List::create(),
    VectorXd downdate_weights = VectorXd::Zero(0), // weights and Xs are considered consecutive. chol_Vi_R will be updated from previous test
    VectorXd inv_prior_X = VectorXd::Zero(0)
) {
  
  int n = Y.rows();
  int m = 1;
  int b_cov = X_cov.cols();
  int b = b_cov + 1;
  int p = X_indices.rows();
  
  bool denseR = Rf_isMatrix(chol_Vi_R_);
  MatrixXd chol_Vi_R_dense;
  SpMat chol_Vi_R_sparse;
  if(denseR){
    chol_Vi_R_dense = as<MatrixXd>(chol_Vi_R_);
  } else {
    chol_Vi_R_sparse = as<SpMat>(chol_Vi_R_);
  }
  
  if(X_cov.rows() != n) stop("Wrong dimenions of X_cov");
  if(X.rows() != n) stop("Wrong dimensions of X");
  if(inv_prior_X.size() < b) {
    if(inv_prior_X.size() == 0) {
      inv_prior_X = VectorXd::Zero(b);
    } else {
      stop("Wrong length of inv_prior_X");
    }
  }
  if(downdate_Xs.size() > 0) {
    if(downdate_Xs.size() < X_indices.col(0).maxCoeff()) stop("Wrong length of downdate_Xs");
  }
  if(X_indices.col(0).maxCoeff() > Y.cols()) stop("X_indices refers to non-existant columns of Y");
  if(X_indices.col(1).maxCoeff() > X.cols()) stop("X_indices refers to non-existant columns of X");
  
  std::vector<SS_result> results;
  
  MatrixXd YX_cov(n,m+b_cov);
  YX_cov.rightCols(b_cov) = X_cov;
  double V_log_det;
  if(denseR) {
    V_log_det = 2*chol_Vi_R_dense.diagonal().array().log().sum();
  } else{
    V_log_det = 2*chol_Vi_R_sparse.diagonal().array().log().sum();
  }
  
  MatrixXd X_cov_std;
  MatrixXd Y_std;
  MatrixXd X_std;
  if(downdate_Xs.size() == 0) {
    if(denseR) {
      Y_std = chol_Vi_R_dense.transpose().triangularView<Lower>().solve(Y);
      X_cov_std = chol_Vi_R_dense.transpose().triangularView<Lower>().solve(X_cov);
      X_std = chol_Vi_R_dense.transpose().triangularView<Lower>().solve(X);
    } else{
      Y_std = chol_Vi_R_sparse.transpose().triangularView<Lower>().solve(Y);
      X_cov_std = chol_Vi_R_sparse.transpose().triangularView<Lower>().solve(X_cov);
      X_std = chol_Vi_R_sparse.transpose().triangularView<Lower>().solve(X);
    }
  } else{
    if(!denseR){
      denseR = true;
      chol_Vi_R_dense = chol_Vi_R_sparse.toDense();
    }
  }
  
  MatrixXd Yi_std, Xi, Xi_std;
  for(int i = 0; i < p; i++) {
    MatrixXd Xi = X.block(0,X_indices.coeff(i,1)-1,n,X_indices.coeff(i,2)); // Note: change to 0-index
    
    int index_i = X_indices.coeff(i,0)-1;
    
    // rotate variables
    if(downdate_Xs.size() > 0 && denseR) {
      // Optionally, downdate chol_Vi_R.
      Rcpp::List downdate_i = as<Rcpp::List>(downdate_Xs[index_i]);  // List with downdate Info
      Rcpp::List downdate_Xis = as<Rcpp::List>(downdate_i["downdate_Xi"]); // List of downdate Matrices (ie for different random effects)
      VectorXd downdate_signs = as<VectorXd>(downdate_i["signs_i"]); // each downdate matrix has same number of columns, and the corrresponding columns are added or removed in each matrix
      int n_downdate = downdate_Xis.length();
      if(n_downdate != downdate_weights.size()) stop("Wrong length of downdate_weights");
      for(int j = 0; j < n_downdate; j++) {
        MatrixXd dXi = as<MatrixXd>(downdate_Xis[j]);
        if(dXi.cols() > 0) {
          if(dXi.rows() != n) stop("Wrong number of rows in downdate_Xs");
          chol_update_R_inplace(chol_Vi_R_dense,dXi,downdate_signs * downdate_weights(j));
        }
      }
      V_log_det = 2*chol_Vi_R_dense.diagonal().array().log().sum();
      
      MatrixXd YX_cov_std;
      YX_cov.leftCols(m) = Y.col(index_i);
      if(denseR) {
        Yi_std = chol_Vi_R_dense.transpose().triangularView<Lower>().solve(Y.col(index_i));
        X_cov_std = chol_Vi_R_dense.transpose().triangularView<Lower>().solve(X_cov);
        Xi_std = chol_Vi_R_dense.transpose().triangularView<Lower>().solve(Xi);
      } else {
        Yi_std = chol_Vi_R_sparse.transpose().triangularView<Lower>().solve(Y.col(index_i));
        X_cov_std = chol_Vi_R_sparse.transpose().triangularView<Lower>().solve(X_cov);
        Xi_std = chol_Vi_R_sparse.transpose().triangularView<Lower>().solve(Xi);
      }
    } else{
      Yi_std = Y_std.col(index_i);
      Xi_std = X_std.block(0,X_indices.coeff(i,1)-1,n,X_indices.coeff(i,2)); // Note: change to 0-index
    }
    
    
    MatrixXd A = X_cov_std.transpose() * X_cov_std;
    A.diagonal() += inv_prior_X.head(b_cov);
    Eigen::LLT<MatrixXd> llt_of_A(A);
    MatrixXd L_A = llt_of_A.matrixL();
    
    if(Xi_std.cols() == 0) {
      // null model
      SS_result result_i = GridLMM_SS(Yi_std,X_cov_std,L_A,V_log_det);
      result_i.ML = calc_ML_c(result_i,n);
      result_i.REML = calc_REML_c(result_i,X_cov);
      results.push_back(result_i);
    } else {
      // Iterate over cis markers
      MatrixXd Xf_std(n,b);
      Xf_std.leftCols(b_cov) = X_cov_std;
      for(int j = 0; j < Xi_std.cols(); j++) {
        MatrixXd Xij_std = Xi_std.col(j);
        
        Xf_std.rightCols(1) = Xij_std;
        
        // calculate L
        MatrixXd B = X_cov_std.transpose() * Xij_std;
        MatrixXd D = Xij_std.transpose() * Xij_std;
        D.diagonal() += inv_prior_X.tail(1);
        MatrixXd L = block_cholesky(L_A, B, D);
        
        // Then solve the model
        SS_result result_i = GridLMM_SS(Yi_std,Xf_std,L,V_log_det);
        
        result_i.ML = calc_ML_c(result_i,n);
        MatrixXd Xf(n,b); 
        Xf << X_cov,Xi.col(j);
        result_i.REML = calc_REML_c(result_i,Xf);
        
        results.push_back(result_i);
      }
    }
  }
  return(collect_SS_results(results));
}
