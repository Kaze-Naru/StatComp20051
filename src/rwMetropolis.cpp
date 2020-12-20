#include <Rcpp.h>
using namespace Rcpp;

//' @title  A random walk Metropolis sampler using Rcpp
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp
//' @param sigma the variance
//' @param x0 the initial value
//' @param N the length of the sample
//' @return a random samples \code{n}
//' @examples
//' \dontrun{
//' rw -> rwMetropolis(2,25,2000)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwMetropolis (double sigma, double x0, int N) {
  NumericVector x(N+1);
  int k = 0;
  x[0] = x0; 
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (exp(-abs(y[0])) / exp(-abs(x[i-1])))){
      x[i] = y[0];
    }
    else { 
      x[i] = x[i-1];
      k=k+1;
    }
  }
  x[N] = k;
  return(x);
} 