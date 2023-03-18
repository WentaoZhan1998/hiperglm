#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
  Eigen::VectorXd QR_solve_rcpp(Eigen::MatrixXd X, Eigen::VectorXd y) {
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(X);
    Eigen::VectorXd beta(qr.solve(y));
    return(beta);
  }
