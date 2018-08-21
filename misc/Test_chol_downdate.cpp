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
  if(L.cols() != X.rows()) stop("dimensions of L and X are not compatible");
  return(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)));
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
      MatrixXd L = block_cholesky2(L_A, B, D);

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
  Rcpp::List proximal_markers
) {
  // Each row of proximal matrix is a vector of 0/1 listing which markers are proximal for a test.
  // X_list is a list of matrices (#=b_x). Each test pulls 1 column from each matrix
  // For each test, we make a make a (n x b_x*l) matrix of of the proximal vectors for the b_x sets of l vectors.
  // The l vectors are only those that differ from the previous test
  // downdate_signs is a list of vectors of +1/-1 for if the corresponding column of the downdate_Xs matrix is to be added or subtracted
                  
  // Note: coming from R, all indexes are 1-based                      
  
  int p = proximal_markers.length();
  
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
  for(int i = 0; i < p; i++){
    // find which markers differ from previous test
    // -1:drop for this test, +1: recover from previous test
    IntegerVector proximal_markers_i = as<IntegerVector>(proximal_markers[i]);
    if(max(proximal_markers_i) > X_list[0].cols()) stop("some proximal markers not present in X_list");
    
    IntegerVector add_markers = setdiff(proximal_markers_i,proximal_markers_0);
    IntegerVector drop_markers = setdiff(proximal_markers_0,proximal_markers_i);
    
    VectorXd downdate_signs_i(add_markers.length() + drop_markers.length());
    downdate_signs_i << VectorXd::Constant(add_markers.length(),-1), VectorXd::Constant(drop_markers.length(),1);
    
    Rcpp::List dXi(n_RE);  // List of downdate matrices for test i (one for each element of X_list)
    Rcpp::List dsi(n_RE);
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
                                        Named("signs_i") = downdate_signs_i);
    
    
    proximal_markers_0 = proximal_markers_i;
  }
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
  int p = downdate_Xs.size();
  
  if(X_cov.rows() != n) stop("Wrong dimenions of X_cov");
  if(inv_prior_X.size() < b) stop("Wrong length of inv_prior_X");
  if(chol_Vi_R.cols() != n || chol_Vi_R.rows() != n) stop("Wrong dimenions of chol_Vi_R");
  if(X_indices.size() != downdate_Xs.length()) stop("Wrong length of X_indices");
  
  std::vector<SS_result> results;
  
  // Process X_list
  std::vector<Map<MatrixXd> > X_list;
  for(int j = 0; j < b_x; j++) {
    X_list.push_back(as<Map<MatrixXd> >(X_list_[j]));
    if(X_list[j].cols() != X_list[0].cols()) stop("Different numbers of columns in X_list matrices");
    if(X_list[j].rows() != n) stop("Wrong number of rows in X_list matrices");
  }
  
  MatrixXd YXf(n,m+b);
  YXf.leftCols(m) = Y;
  YXf.block(0,m,n,b_cov) = X_cov;
  for(int i = 0; i < p; i++) {

    // update chol_Vi_R and V_log_det

    Rcpp::List downdate_i = as<Rcpp::List>(downdate_Xs[i]);  // List with downdate Info
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
      int index_i = X_indices[i]-1;
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

    results.push_back(result_i);
  }

  return(collect_SS_results(results));
}

Rcpp::List GridLMM_SSb(
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
    VectorXd V_star_L = vectorize_chol2(L_A);
    VectorXd V_star_inv_log_det(1);
    V_star_inv_log_det(0) = 2*L_A.diagonal().array().log().sum();
    MatrixXd b_hat = solve_cholesky2(L_A,X_cov_std.transpose() * Y_std);
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
    MatrixXd L = block_cholesky2(L_A, B, D);
    
    V_star_L.col(i) = vectorize_chol2(L);
    V_star_inv_log_det(i) = 2*L.diagonal().array().log().sum();
    
    MatrixXd X_std(n,b);
    X_std << X_cov_std,Xi_std;
    MatrixXd b_hat = solve_cholesky2(L,X_std.transpose() * Y_std);
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
Rcpp::List GridLMM_SS_dense_c3(
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
      // Xi.col(j) = inv_chol_Vi_transpose.triangularView<Lower>() * Xi_full.col(X_indices[j]-1);
    }
    X_stds.push_back(inv_chol_Vi_transpose.triangularView<Lower>() * Xi);
    // X_stds.push_back(Xi);
    // X_stds.push_back(inv_chol_Vi_transpose * Xi);
  }
  
  MatrixXd X_cov_std = inv_chol_Vi_transpose.triangularView<Lower>() * X_cov;
  
  return(GridLMM_SSb(Y_std,X_cov_std,X_stds,inv_prior_X,V_log_det));
}