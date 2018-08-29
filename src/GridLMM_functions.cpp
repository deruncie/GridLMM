#include <math.h>
#include "GridLMM_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// [[Rcpp::export()]]
Rcpp::List svd_c(Map<MatrixXd> X) {
  // Replacement of R's svd function with Eigen's version
  Eigen::BDCSVD<MatrixXd> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
  return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
                            Named("d") = bcdsolve.singularValues(),
                            Named("v") = bcdsolve.matrixV()));
}


// [[Rcpp::export()]]
MatrixXd chol_c(Map<MatrixXd> X){
  // Replacement of R's chol function with Eigen's version
  Eigen::LLT<MatrixXd> llt_of_X(X);
  MatrixXd U = llt_of_X.matrixU();
  return(U);
}

// [[Rcpp::export()]]
Rcpp::List premultiply_list_of_matrices(SEXP Qt_, Rcpp::List X_list){
  int p = X_list.length();
  Rcpp::List X_list_new;
  
  bool Qt_dense = Rf_isMatrix(Qt_);
  MatrixXd Qt_matrix;
  if(Qt_dense) {
    Qt_matrix = as<MatrixXd>(Qt_);
  } 
  
  for(int i = 0; i < p; i++){
    MatrixXd Xi = Rcpp::as<Map<MatrixXd> >(X_list[i]);
    if(Qt_dense) {
      X_list_new.push_back(Qt_matrix * Xi);
    } else{
      MSpMat Qt_SpMat = as<MSpMat>(Qt_);
      X_list_new.push_back(Qt_SpMat * Xi);
    }
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


MatrixXd block_cholesky(MatrixXd L_A, // Cholesky L matrix for A matrix (lower-left block)
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

void chol_update_R_inplace(MatrixXd &R, MatrixXd X, VectorXd weights) {
  int n = R.rows();
  if(X.rows() != n) stop("Wrong dimension of X for downdating R");
  if(weights.size() != X.cols()) stop("wrong length of weights for downdating R");
  for(int i = 0; i < X.cols(); i++){
    VectorXd Xi = X.col(i);
    for(int k = 0; k < n; k++){
      double R_kk = R.coeffRef(k,k);
      double x_k = Xi[k];
      double r = sqrt(R_kk*R_kk + weights[i] * x_k*x_k);
      double c = r / R_kk;
      double s = x_k / R_kk;
      R.coeffRef(k,k) = r;
      if(k < (n-1)) {
        R.block(k,k+1,1,n-k-1) = (R.block(k,k+1,1,n-k-1) + weights[i] * s * Xi.tail(n-k-1).transpose())/c;
        Xi.tail(n-k-1) = c*Xi.tail(n-k-1) - s*R.block(k,k+1,1,n-k-1).transpose();
      }
    }
  }
}


// [[Rcpp::export()]]
MatrixXd chol_update_L(MatrixXd L, MatrixXd X, VectorXd weights) {
  MatrixXd R = L.transpose();
  chol_update_R_inplace(R,X,weights);
  return(R.transpose());
}

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

// [[Rcpp::export()]]
MatrixXd crossprod_cholR(Map<MatrixXd> chol_R, Map<MatrixXd> X){
  // t(chol_R) %*% chol_R = K
  // returns t(chol_R) %*% X
  return(chol_R.transpose().triangularView<Lower>() * X);
}


ArrayXd calc_ML_c(SS_result SS, int n){
  double c = n*std::log(n/(2.0*M_PI)) - n - SS.V_log_det;
  ArrayXd ML = -n/2.0*SS.RSSs.array().log() + c/2.0;
  return(ML);
}

ArrayXd calc_REML_c(SS_result SS, MatrixXd &X){
  int n = X.rows();
  double b = X.cols();
  // int m = SS.RSSs.size();
  ArrayXd ML = SS.ML;
  if(SS.ML.size() == 0) {
    ML = calc_ML_c(SS,n);
  }
  
  MatrixXd XtX = X.transpose() * X;
  Eigen::LLT<MatrixXd> llt_of_XtX(XtX);
  MatrixXd L = llt_of_XtX.matrixL();
  double log_det_X = 2*L.diagonal().array().log().sum();
  
  ArrayXd REML = ML + 0.5 * (b*(2.0*M_PI/(n-b)*SS.RSSs.array()).log() + log_det_X - SS.V_star_inv_log_det);
  return(REML);
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

SS_result GridLMM_SS(
    MatrixXd &Y_std,
    MatrixXd &X_std,
    MatrixXd &L,
    double V_log_det
) {
  SS_result result;
  
  result.V_log_det = V_log_det;
  result.V_star_L = vectorize_chol(L);
  result.V_star_inv_log_det = 2*L.diagonal().array().log().sum();
  
  MatrixXd b_hat = solve_cholesky(L,X_std.transpose() * Y_std);
  MatrixXd eta_std = Y_std - X_std * b_hat;
  
  result.beta_hats = b_hat;
  result.RSSs = eta_std.cwiseProduct(eta_std).colwise().sum();
  
  result.ML = ArrayXd::Constant(Y_std.cols(),0);
  result.REML = ArrayXd::Constant(Y_std.cols(),0);
  
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
  
  MatrixXd MLs(m,p);
  MatrixXd REMLs(m,p);
  
  for(int i = 0; i < p; i++){
    SS_result results_i = results[i];
    beta_hats.col(i) = Map<VectorXd>(results_i.beta_hats.data(),results_i.beta_hats.size());
    RSSs.col(i) = results_i.RSSs;
    V_star_L.col(i) = results_i.V_star_L;
    V_star_inv_log_det(i) = results_i.V_star_inv_log_det;
    V_log_det(i) = results_i.V_log_det;
    
    MLs.col(i) = results_i.ML;
    REMLs.col(i) = results_i.REML;
  }
  
  return(Rcpp::List::create(Named("beta_hats") = beta_hats,
                            Named("RSSs") = RSSs,
                            Named("V_log_dets") = V_log_det,
                            Named("V_star_inv_log_det") = V_star_inv_log_det,
                            Named("V_star_L") = V_star_L
                            ,Named("MLs") = MLs
                            ,Named("REMLs") = REMLs
                            )
           );
}


// [[Rcpp::export()]]
Rcpp::List GridLMM_SS_matrix(
    Map<MatrixXd> Y,
    SEXP chol_Vi_R_,
    Map<MatrixXd> X_cov,
    Rcpp::List X_list_,
    ArrayXi X_indices,
    VectorXd inv_prior_X
) {
  
  bool denseR = Rf_isMatrix(chol_Vi_R_);
  
  int n = Y.rows();
  int b_cov = X_cov.cols();
  int b_x = X_list_.size();
  int b = b_cov + b_x;
  int p = X_indices.size();
  
  if(X_cov.rows() != n) stop("Wrong dimenions of X_cov");
  if(inv_prior_X.size() < b) stop("Wrong length of inv_prior_X");
  
  std::vector<SS_result> results;
  
  // Process X_list
  std::vector<Map<MatrixXd> > X_list;
  for(int j = 0; j < b_x; j++) {
    X_list.push_back(as<Map<MatrixXd> >(X_list_[j]));
    if(X_list[j].cols() != X_list[0].cols()) stop("Different numbers of columns in X_list matrices");
    if(X_list[j].rows() != n) stop("Wrong number of rows in X_list matrices");
  }
  
  if(b_x > 0) {
    if(X_indices.maxCoeff() > X_list[0].cols()) stop("X_indices reference missing columns of X_list");
    if(X_indices.minCoeff() < 1) stop("X_indices must be >= 1");
  }
  
  MatrixXd Y_std, X_cov_std, X_std;
  double V_log_det = 0;
  if(denseR) {
    Map<MatrixXd> chol_Vi_R = as<Map<MatrixXd> >(chol_Vi_R_);
    if(chol_Vi_R.cols() != n || chol_Vi_R.rows() != n) stop("Wrong dimenions of chol_Vi_R");
    V_log_det = 2*chol_Vi_R.diagonal().array().log().sum();
    Y_std = chol_Vi_R.transpose().triangularView<Lower>().solve(Y);
    X_cov_std = chol_Vi_R.transpose().triangularView<Lower>().solve(X_cov);
  } else{
    MSpMat chol_Vi_R = as<MSpMat>(chol_Vi_R_);
    if(chol_Vi_R.cols() != n || chol_Vi_R.rows() != n) stop("Wrong dimenions of chol_Vi_R");
    for(int i = 0; i < n; i++) V_log_det += 2*log(chol_Vi_R.coeff(i,i));
    Y_std = chol_Vi_R.transpose().triangularView<Lower>().solve(Y);
    X_cov_std = chol_Vi_R.transpose().triangularView<Lower>().solve(X_cov);
  }
  
  MatrixXd A = X_cov_std.transpose() * X_cov_std;
  A.diagonal() += inv_prior_X.head(b_cov);
  Eigen::LLT<MatrixXd> llt_of_A(A);
  MatrixXd L_A = llt_of_A.matrixL();
  
  if(b_x == 0) {
    // If there are no tests to do, just run the null model
    MatrixXd Xf_std = X_cov_std;
    SS_result result_1 = GridLMM_SS(Y_std,Xf_std,L_A,V_log_det);
    results.push_back(result_1);
  } else{
    // process X for active_X_columns
    MatrixXd X_std(n,p*b_x);
    for(int i = 0; i < p; i++) {
      // reform X_list into a wide matrix of n x b_x design matrices
      int index_i = X_indices[i]-1;
      for(int j = 0; j < b_x; j++) {
        X_std.col(i*b_x + j) = X_list[j].col(index_i);
      }
    }
    // rotate these design matrixes by chol_Vi_R
    if(denseR) {
      Map<MatrixXd> chol_Vi_R = as<Map<MatrixXd> >(chol_Vi_R_);
      X_std = chol_Vi_R.transpose().triangularView<Lower>().solve(X_std);
    } else{
      MSpMat chol_Vi_R = as<MSpMat>(chol_Vi_R_);
      X_std = chol_Vi_R.transpose().triangularView<Lower>().solve(X_std);
    }
    
    // run tests
    for(int i = 0; i < p; i++) {
      MatrixXd Xi_std = X_std.block(0,i*b_x,n,b_x);
      
      MatrixXd Xf_std(n,b);
      Xf_std << X_cov_std,Xi_std;
      
      // calculate L
      MatrixXd B = X_cov_std.transpose() * Xi_std;
      MatrixXd D = Xi_std.transpose() * Xi_std;
      D.diagonal() += inv_prior_X.tail(b_x);
      MatrixXd L = block_cholesky(L_A, B, D);
      
      // Then solve the model
      SS_result result_i = GridLMM_SS(Y_std,Xf_std,L,V_log_det);
      
      results.push_back(result_i);
    }
  }
  
  return(collect_SS_results(results));
}

// [[Rcpp::export()]]
Rcpp::List build_downdate_Xs(
    IntegerVector RE_index,
    Rcpp::List X_list_,
    Rcpp::List proximal_markers,
    IntegerVector active_X_list
) {
  // Each row of proximal matrix is a vector of 0/1 listing which markers are proximal for a test.
  // X_list is a list of matrices (#=b_x). Each test pulls 1 column from each matrix
  // For each test, we make a make a (n x b_x*l) matrix of of the proximal vectors for the b_x sets of l vectors.
  // The l vectors are only those that differ from the previous test
  // downdate_signs is a list of vectors of +1/-1 for if the corresponding column of the downdate_Xs matrix is to be added or subtracted
  // downdate_Xs inherits proximal_markers names
  
  // Note: coming from R, all indexes are 1-based
  
  int p = proximal_markers.length();
  int p_active = active_X_list.size();
  if(max(active_X_list) > p) stop("active_X_list refers to non-existant proximal_markers list");
  
  Rcpp::List downdate_Xs(p);
  
  // Process X_list
  RE_index = RE_index - 1;  // change to 0-based
  int n_RE = RE_index.length();
  if(max(RE_index) > X_list_.length()) stop("RE_index refers to non-existant X-matrix");
  std::vector<Map<MatrixXd> > X_list;
  for(int j = 0; j < n_RE; j++) {
    if(RE_index(j) == -1) {
      MatrixXd null = MatrixXd::Zero(0,0);
      X_list.push_back(Map<MatrixXd>(null.data(),0,0));
    } else{
      X_list.push_back(as<Map<MatrixXd> >(X_list_[RE_index(j)]));
    }
  }
  int n = X_list[0].rows();
  
  IntegerVector proximal_markers_0(0);
  for(int i_active = 0; i_active < p_active; i_active++){
    int i = active_X_list(i_active)-1;  // i indexes proximal_markers (full marker order). i_active indexes active_X_list
    // find which markers differ from previous test
    // -1:drop for this test, +1: recover from previous test
    IntegerVector proximal_markers_i(0);
    if(Rf_isInteger(proximal_markers[i])) {
      proximal_markers_i = as<IntegerVector>(proximal_markers[i]);
    }
    if(max(proximal_markers_i) > X_list[0].cols()) stop("some proximal markers not present in X_list");

    IntegerVector add_markers = setdiff(proximal_markers_i,proximal_markers_0);
    IntegerVector drop_markers = setdiff(proximal_markers_0,proximal_markers_i);

    VectorXd downdate_signs_i(add_markers.length() + drop_markers.length());
    downdate_signs_i << VectorXd::Constant(add_markers.length(),-1), VectorXd::Constant(drop_markers.length(),1);

    Rcpp::List dXi(n_RE);  // List of downdate matrices for test i (one for each element of X_list)

    for(int j = 0; j < n_RE; j++) {
      if(X_list[j].cols() > 0) {
        MatrixXd dXij(n,add_markers.length() + drop_markers.length());
        for(int k = 0; k < add_markers.length(); k++) {
          dXij.col(k) = X_list[j].col(add_markers(k)-1);
        }
        for(int k = 0; k < drop_markers.length(); k++) {
          dXij.col(add_markers.length() + k) = X_list[j].col(drop_markers(k)-1);
        }
        dXi[j] = dXij;
      } else{
        dXi[j] = MatrixXd::Zero(n,0);
      }
    }
    downdate_Xs[i] = Rcpp::List::create(Named("downdate_Xi") = dXi,
                                        Named("signs_i") = downdate_signs_i,
                                        Named("add_markers") = add_markers,
                                        Named("drop_markers") = drop_markers);
    
    
    proximal_markers_0 = proximal_markers_i;
  }
  downdate_Xs.names() = proximal_markers.names();
  return(downdate_Xs);
}

// [[Rcpp::export()]]
Rcpp::List GridLMM_SS_downdate_matrix(
    Map<MatrixXd> Y,
    MatrixXd chol_Vi_R, //don't want this modified
    Map<MatrixXd> X_cov,
    Rcpp::List X_list_,
    ArrayXi X_indices,
    Rcpp::List downdate_Xs,
    VectorXd downdate_weights, // weights and Xs are considered consecutive. chol_Vi_R will be updated from previous test
    VectorXd inv_prior_X
) {
  
  int n = Y.rows();
  int m = Y.cols();
  int b_cov = X_cov.cols();
  int b_x = X_list_.length();
  int b = b_cov + b_x;
  int p = X_indices.size();
  
  if(X_cov.rows() != n) stop("Wrong dimenions of X_cov");
  if(inv_prior_X.size() < b) stop("Wrong length of inv_prior_X");
  if(chol_Vi_R.cols() != n || chol_Vi_R.rows() != n) stop("Wrong dimenions of chol_Vi_R");
  if(X_indices.maxCoeff() > downdate_Xs.length()) stop("Bad indices in X_indices");
  
  std::vector<SS_result> results;
  
  // Process X_list
  std::vector<Map<MatrixXd> > X_list;
  for(int j = 0; j < b_x; j++) {
    X_list.push_back(as<Map<MatrixXd> >(X_list_[j]));
    if(X_list[j].cols() != downdate_Xs.length()) stop("Number of columns in X_list matrices not equal to length of downdate_Xs");
    if(X_list[j].rows() != n) stop("Wrong number of rows in X_list matrices");
  }
  
  MatrixXd YXf(n,m+b);
  YXf.leftCols(m) = Y;
  YXf.block(0,m,n,b_cov) = X_cov;
  for(int i = 0; i < p; i++) {
    int index_i = X_indices[i]-1;
    
    // update chol_Vi_R and V_log_det
    
    Rcpp::List downdate_i = as<Rcpp::List>(downdate_Xs[index_i]);  // List with downdate Info
    Rcpp::List downdate_Xis = as<Rcpp::List>(downdate_i["downdate_Xi"]); // List of downdate Matrices (ie for different random effects)
    VectorXd downdate_signs = as<VectorXd>(downdate_i["signs_i"]); // each downdate matrix has same number of columns, and the corrresponding columns are added or removed in each matrix
    int n_downdate = downdate_Xis.length();
    if(n_downdate != downdate_weights.size()) stop("Wrong length of downdate_weights");
    for(int j = 0; j < n_downdate; j++) {
      MatrixXd dXi = as<MatrixXd>(downdate_Xis[j]);
      if(dXi.cols() > 0) {
        if(dXi.rows() != n) stop("Wrong number of rows in downdate_Xs");
        chol_update_R_inplace(chol_Vi_R,dXi,downdate_signs * downdate_weights(j));
      }
    }
    double V_log_det = 2*chol_Vi_R.diagonal().array().log().sum();
    
    // re-form column of X into a n x b_x design matrix
    if(b_x > 0) {
      MatrixXd Xi(n,b_x);
      for(int j = 0; j < b_x; j++) {
        Xi.col(j) = X_list[j].col(index_i);
      }
      YXf.rightCols(b_x) = Xi;
    }
    
    // rotate YXf matrix, extract Y_std and Xf_std;
    MatrixXd YXf_std = chol_Vi_R.transpose().triangularView<Lower>().solve(YXf);
    MatrixXd Y_std = YXf_std.leftCols(m);
    MatrixXd Xf_std = YXf_std.rightCols(b);
    
    // calculate L
    MatrixXd XtX = Xf_std.transpose() * Xf_std;
    XtX.diagonal() += inv_prior_X;
    Eigen::LLT<MatrixXd> llt_of_XtX(XtX);
    MatrixXd L = llt_of_XtX.matrixL();
    
    // Then solve the model
    SS_result result_i = GridLMM_SS(Y_std,Xf_std,L,V_log_det);
    
    
    // MatrixXd Xi = YXf.rightCols(b);
    // result_i.ML = calc_ML_c(result_i,Y_std.rows());
    // result_i.REML = calc_REML_c(result_i,Xi);
    
    results.push_back(result_i);
  }
  
  return(collect_SS_results(results));
}
