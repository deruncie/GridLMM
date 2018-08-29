#include "GridLMM_types.h"
#include "brent.hpp"
using namespace Rcpp;

// This is more efficient for updating
void chol_update_L_inplace(MatrixXd &L, MatrixXd X, VectorXd weights) {
  int n = L.rows();
  if(X.rows() != n) stop("Wrong dimension of X for updating L");
  if(weights.size() != X.cols()) stop("wrong length of weights for updating L");
  for(int i = 0; i < X.cols(); i++){
    VectorXd Xi = X.col(i);
    double weights_i = weights[i];
    for(int k = 0; k < n; k++){
      double x_k = Xi.coeffRef(k);
      if(x_k != 0) {
        double L_kk = L.coeffRef(k,k);
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
        }
      }
    }
  }
}


double fit_model_tau(SS_result &SS,double tau, MatrixXd &YXf, MatrixXd &G, 
                     MatrixXd chol_Vi_L, VectorXd &inv_prior_X, int m, int b) {
  
  int k = G.cols();
  int n = YXf.rows();
  
  VectorXd weights = VectorXd::Constant(k,tau/(1-tau)*1/k);
  // Rcout << 1 << " ";
  chol_update_L_inplace(chol_Vi_L,G,weights);
  // Rcout << 2 << " ";
  chol_Vi_L *= std::sqrt(1-tau);
  double V_log_det = 2*chol_Vi_L.diagonal().array().log().sum();
  // Rcout << V_log_det << " " << tau  << " ";
  
  // rotate YXf matrix, extract Y_std and Xf_std;
  // Rcout << 3 << " ";
  MatrixXd YXf_std = chol_Vi_L.triangularView<Lower>().solve(YXf);
  MatrixXd Y_std = YXf_std.leftCols(m);
  MatrixXd Xf_std = YXf_std.rightCols(b);
  
  // calculate L
  // Rcout << 4 << " ";
  MatrixXd XtX = Xf_std.transpose() * Xf_std;
  XtX.diagonal() += inv_prior_X;
  Eigen::LLT<MatrixXd> llt_of_XtX(XtX);
  MatrixXd L = llt_of_XtX.matrixL();
  
  // Then solve the model
  
  // Rcout << 5 << " ";
  SS = GridLMM_SS(Y_std,Xf_std,L,V_log_det);
  // Rcout << SS.RSSs << std::endl;
  
  MatrixXd X = YXf.rightCols(b);
  // Rcout << 6 << " ";
  double REML = calc_REML_c(SS,X).coeff(0);
  // Rcout << "done\n";
  return(REML);
}


class optimizerTau : public brent::func_base
{
  private:
    MatrixXd &YXf;
    MatrixXd &G;
    MatrixXd &chol_Vi_L;
    SS_result &SS;
    VectorXd &inv_prior_X;
    int m,b;
    int n_calls;
  public:
    optimizerTau (MatrixXd &YXf_, MatrixXd &G_, MatrixXd &chol_Vi_L_, SS_result &SS_,
                 VectorXd &inv_prior_X_, int m_, int b_) : 
      YXf(YXf_), G(G_), chol_Vi_L(chol_Vi_L_), SS(SS_),
      inv_prior_X(inv_prior_X_), m(m_), b(b_)
    {
      n_calls = 0;
    }
    int get_n_calls() {
      return(n_calls);  
    }
  double operator() (double tau) {
    if(m != 1) {
      stop("fitting only works with m=1");
    }
    n_calls++;
    double REML = fit_model_tau(SS, tau, YXf, G, chol_Vi_L, inv_prior_X, m, b);
    return(-REML);
  }
};


// [[Rcpp::export()]]
Rcpp::List GridLMM_test_setTest(
  MatrixXd Y,
  MatrixXd chol_Vi_R,
  MatrixXd G,
  MatrixXd X,
  double tolerance = 0.01
) {
  int n = Y.rows();
  int m = Y.cols();
  int b = X.cols();
  MatrixXd YXf(n,m+b);
  YXf << Y,X;
  
  VectorXd inv_prior_X = VectorXd::Zero(b);
  
  // optimize tau
  SS_result result_i;
  MatrixXd chol_Vi_L = chol_Vi_R.transpose();
  optimizerTau get_tau(YXf, G, chol_Vi_L, result_i,
                       inv_prior_X, m,b);
  
  double tau;
  brent::local_min(0.0, 0.99, tolerance, get_tau, tau); // maximum of tau is 0.9
  
  double REML = fit_model_tau(result_i,tau,YXf, G, chol_Vi_L, inv_prior_X, m,b);
  
  
  std::vector<SS_result> results;
  results.push_back(result_i);
  Rcpp::List SS_results = collect_SS_results(results);
  SS_results["REMLs"] = REML;
  SS_results["taus"] = tau;
  SS_results["n_calls"] = get_tau.get_n_calls();
  return(SS_results);
}

// [[Rcpp::export()]]
Rcpp::List GridLMM_setTest_downdate_matrix(
    Map<MatrixXd> Y,
    MatrixXd chol_Vi_R, //don't want this modified
    Map<MatrixXd> X_cov,
    Map<MatrixXd> X,
    ArrayXXi X_indices, // p x 3 array. 1st column is test index (1-based, corresponds to downdate_Xs). columns 2 and 3 are start and end (1-based) of set
    Rcpp::List downdate_Xs,
    VectorXd downdate_weights, // weights and Xs are considered consecutive. chol_Vi_R will be updated from previous test
    VectorXd inv_prior_X,
    double tolerance = 0.01
) {
  
  int n = Y.rows();
  int m = Y.cols();
  int b = X_cov.cols();
  int p = X_indices.rows();
  
  if(X_cov.rows() != n) stop("Wrong dimenions of X_cov");
  if(inv_prior_X.size() < b) stop("Wrong length of inv_prior_X");
  if(chol_Vi_R.cols() != n || chol_Vi_R.rows() != n) stop("Wrong dimenions of chol_Vi_R");
  if(X_indices.col(0).maxCoeff() > downdate_Xs.length()) stop("Bad indices in X_indices");
  
  std::vector<SS_result> results;
  VectorXd taus(p);
  VectorXd REMLs(p);
  
  MatrixXd YXf(n,m+b);
  YXf.leftCols(m) = Y;
  YXf.rightCols(b) = X_cov;
  for(int i = 0; i < p; i++) {
    int downdate_index = X_indices.coeff(i,0) - 1;
    int set_start = X_indices.coeff(i,1) - 1;
    int set_length = X_indices.coeff(i,2) - set_start;
    
    // update chol_Vi_R based on downdate_Xs
    
    Rcpp::List downdate_i = as<Rcpp::List>(downdate_Xs[downdate_index]);  // List with downdate Info
    Rcpp::List downdate_Xis = as<Rcpp::List>(downdate_i["downdate_Xi"]); // List of downdate Matrices (ie for different random effects)
    VectorXd downdate_signs = as<VectorXd>(downdate_i["signs_i"]); // each downdate matrix has same number of columns, and the corrresponding columns are added or removed in each matrix
    int n_downdate = downdate_Xis.length();
    if(n_downdate != 1) stop("Only one downdate matrix allowed");
    for(int j = 0; j < n_downdate; j++) {
      MatrixXd dXi = as<MatrixXd>(downdate_Xis[j]);
      if(dXi.cols() > 0) {
        if(dXi.rows() != n) stop("Wrong number of rows in downdate_Xs");
        chol_update_R_inplace(chol_Vi_R,dXi,downdate_signs * downdate_weights(j));
      }
    }
    
    double tau;
    double REML;
    MatrixXd G;
    SS_result result_i;
    MatrixXd chol_Vi_L = chol_Vi_R.transpose();
    if(X.cols() > 0) {
      // Extract G matrix
      G = X.block(0,set_start,n,set_length);
      
      // copy chol_Vi_L for optimizing tau
      
      // optimize tau
      optimizerTau get_tau(YXf, G, chol_Vi_L, result_i,
                           inv_prior_X, m,b);
      brent::local_min(0.0, 0.99, tolerance, get_tau, tau); // maximum of tau is 0.99
    } else{
      G = MatrixXd::Zero(n,0);
      tau = 0;
    }
    REML = fit_model_tau(result_i,tau,YXf, G, chol_Vi_L, inv_prior_X, m,b);
    REMLs(i) = REML;
    taus(i) = tau;
    results.push_back(result_i);
  }
  Rcpp::List SS_results = collect_SS_results(results);
  SS_results["REML"] = REMLs;
  SS_results["tau"] = taus;
    
  return(SS_results);
}







double fit_model_tau2(SS_result &SS,double tau, MatrixXd &YXf, VectorXd G_d, 
                     double V_log_det_R, MatrixXd &X, VectorXd &inv_prior_X, int m, int b) {
  
  int k = G_d.size();
  int n = YXf.rows();
  
  // double weight = tau/(1-tau)*1/k;
  
  VectorXd diagV = VectorXd::Constant(n,1.0-tau);
  diagV.head(k) += tau/k * G_d.cwiseProduct(G_d);
  
  double V_log_det = V_log_det_R + diagV.array().log().sum();
  // Rcout << V_log_det << " " << tau  << " ";
  
  // rotate YXf matrix, extract Y_std and Xf_std;
  // Rcout << 3 << " ";
  MatrixXd YXf_std = diagV.cwiseSqrt().cwiseInverse().asDiagonal() * YXf;
  MatrixXd Y_std = YXf_std.leftCols(m);
  MatrixXd Xf_std = YXf_std.rightCols(b);
  
  // calculate L
  // Rcout << 4 << " ";
  MatrixXd XtX = Xf_std.transpose() * Xf_std;
  XtX.diagonal() += inv_prior_X;
  Eigen::LLT<MatrixXd> llt_of_XtX(XtX);
  MatrixXd L = llt_of_XtX.matrixL();
  
  // Then solve the model
  
  // Rcout << 5 << " ";
  SS = GridLMM_SS(Y_std,Xf_std,L,V_log_det);
  // Rcout << SS.RSSs<< std::endl;
  
  // MatrixXd X = YXf.rightCols(b);
  // Rcout << 6 << " ";
  double REML = calc_REML_c(SS,X).coeff(0);
  // Rcout << "done\n";
  return(REML);
}


class optimizerTau2 : public brent::func_base
{
private:
  MatrixXd &YXf;
  VectorXd &G_d;
  double V_log_det_R;
  MatrixXd &X;
  VectorXd &inv_prior_X;
  int m,b;
  SS_result &SS;
  int n_calls;
public:
  optimizerTau2 (MatrixXd &YXf_, VectorXd &G_d_, double V_log_det_R_, MatrixXd &X_,
                VectorXd &inv_prior_X_, int m_, int b_, SS_result &SS_) : 
  YXf(YXf_), G_d(G_d_), V_log_det_R(V_log_det_R_), X(X_),
  inv_prior_X(inv_prior_X_), m(m_), b(b_), SS(SS_)
  {
    n_calls = 0;
  }
  int get_n_calls() {
    return(n_calls);  
  }
  double operator() (double tau) {
    if(m != 1) {
      stop("fitting only works with m=1");
    }
    n_calls++;
    double REML = fit_model_tau2(SS, tau, YXf, G_d, V_log_det_R, X, inv_prior_X, m, b);
    return(-REML);
  }
};


Rcpp::List svd_c2(MatrixXd &X) {
  // Replacement of R's svd function with Eigen's version
  Eigen::BDCSVD<MatrixXd> bcdsolve(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
  return(Rcpp::List::create(Named("u") = bcdsolve.matrixU(),
                            Named("d") = bcdsolve.singularValues(),
                            Named("v") = bcdsolve.matrixV()));
}

// [[Rcpp::export()]]
Rcpp::List GridLMM_test_setTest2(
    MatrixXd Y,
    MatrixXd chol_Vi_R,
    MatrixXd G,
    MatrixXd X,
    double tolerance = 0.01
) {
  int n = Y.rows();
  int m = Y.cols();
  int b = X.cols();
  
  double V_log_det_R = 2*chol_Vi_R.diagonal().array().log().sum();
  
  MatrixXd Rinvt_G = chol_Vi_R.transpose().triangularView<Lower>().solve(G);
  MatrixXd YXf(n,m+b);
  YXf << Y,X;
  MatrixXd YXf_std = chol_Vi_R.transpose().triangularView<Lower>().solve(YXf);
  
  Rcpp::List svd_Rinvt_G = svd_c2(Rinvt_G);
  MatrixXd U = as<Map<MatrixXd> >(svd_Rinvt_G["u"]);
  VectorXd G_d = as<VectorXd>(svd_Rinvt_G["d"]);
  MatrixXd Ut_YXf_std = U.transpose() * YXf_std;
  
  VectorXd inv_prior_X = VectorXd::Zero(b);
  
  // optimize tau
  SS_result result_i;
  MatrixXd chol_Vi_L = chol_Vi_R.transpose();
  optimizerTau2 get_tau(Ut_YXf_std, G_d, V_log_det_R,X,
                       inv_prior_X, m,b,result_i);
  
  double tau;
  brent::local_min(0.0, 0.99, tolerance, get_tau, tau); // maximum of tau is 0.9
  
  double REML = fit_model_tau2(result_i,tau,Ut_YXf_std, G_d, V_log_det_R,X,inv_prior_X, m,b);
  
  
  std::vector<SS_result> results;
  results.push_back(result_i);
  Rcpp::List SS_results = collect_SS_results(results);
  SS_results["REMLs"] = REML;
  SS_results["taus"] = tau;
  SS_results["n_calls"] = get_tau.get_n_calls();
  return(SS_results);
}









// // [[Rcpp::export]]
// double f(double x){
//   myFunctorClass my_f(4.0, 1.0);
//   return my_f(x);
// }
// 
// // [[Rcpp::export]]
// List min_f(){
//   myFunctorClass my_f(4.0, 1.0);
//   double x_star;
//   double f_star = brent::local_min(0.0, 8.0, 0.0001, my_f, x_star);
//   return List::create(Named("minimum") = x_star, Named("objective") = f_star);
// }