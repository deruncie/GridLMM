
#include "GridLMM_types.h"
#include "brent.hpp"
using namespace Rcpp;

class myFunctorClass : public brent::func_base
{
  private:
    double a, b;
  public:
    myFunctorClass (double a_, double b_) : a(a_), b(b_) {}
  double operator() (double x) {
    return (x - a) * (x - a) + b;
  }
};

// [[Rcpp::export]]
double f(double x){
  myFunctorClass my_f(4.0, 1.0);
  return my_f(x);
}

// [[Rcpp::export]]
List min_f(){
  myFunctorClass my_f(4.0, 1.0);
  double x_star;
  double f_star = brent::local_min(0.0, 8.0, 0.0001, my_f, x_star);
  return List::create(Named("minimum") = x_star, Named("objective") = f_star);
}