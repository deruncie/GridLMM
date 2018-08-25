#include <RcppEigen.h>


using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::ArrayXd;                  // variable size vector, double precision
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;

using namespace Rcpp;
using namespace Eigen;


// Helper function and class definitions

struct SS_result{
  MatrixXd beta_hats;
  VectorXd RSSs;
  double V_log_det;
  double V_star_inv_log_det;
  VectorXd V_star_L;
  ArrayXd ML,REML;
};

Rcpp::List svd_c(Map<MatrixXd> X);
MatrixXd chol_c(Map<MatrixXd> X);
Rcpp::List premultiply_list_of_matrices(SEXP Qt_, Rcpp::List X_list);
double log_det_from_LDLt(Eigen::LDLT<MatrixXd> ldlt);
MatrixXd block_cholesky(MatrixXd L_A, MatrixXd B, MatrixXd D);
MatrixXd solve_cholesky(MatrixXd L, MatrixXd X);
VectorXd vectorize_chol(MatrixXd L);
MatrixXd fill_chol_L(VectorXd vec_chol,int p);
void chol_update_R_inplace(MatrixXd &R, MatrixXd X, VectorXd weights);
MatrixXd chol_update_L(MatrixXd L, MatrixXd X, VectorXd weights);
MatrixXd chol_update(MatrixXd L, MatrixXd X, int sign);
MatrixXd chol_dropRows(Map<MatrixXd> L,int start_row, int num_rows);
MatrixXd crossprod_cholR(Map<MatrixXd> chol_R, Map<MatrixXd> X);
MatrixXd F_hats(Map<MatrixXd> beta_hats,Map<MatrixXd> RSSs,Map<MatrixXd> V_star_L,int n,int b,int m);
VectorXd log_det_of_XtX(Map<MatrixXd> X_cov,Rcpp::List X_tests,ArrayXi X_indices);

SS_result GridLMM_SS(MatrixXd &Y_std,MatrixXd &X_std,MatrixXd &L,double V_log_det);
Rcpp::List collect_SS_results(std::vector<SS_result> results);
Rcpp::List build_downdate_Xs(IntegerVector RE_index,Rcpp::List X_list_,Rcpp::List proximal_markers,IntegerVector active_X_list);

ArrayXd calc_ML_c(SS_result SS, int n);
ArrayXd calc_REML_c(SS_result SS, MatrixXd &X);
