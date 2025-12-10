#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// 4/5/2023 debug:overlapping variable name
//
// pooi = pool index
// pind: current individual probability 
// z = pool response

// [[Rcpp::export]]
NumericVector sampletildey_it(NumericVector pind, NumericVector yinit,
                              double Se, double Sp){
  int nsample = yinit.length();
  NumericVector ysample(nsample);
  
  for(int i=0; i<nsample; ++i){
    double pind0 = pind[i]; // probability of individual i
    double z0 = yinit[i]; // result for individual testing
    double p1 = pind0 * pow(Se,z0) * pow(1-Se, 1-z0);
    double z1 = 0; //new subset pool response for individual testing
    double p0 = (1-pind0) * pow( pow(Se, z0)*pow(1-Se, 1-z0), z1) *
      pow( pow(1-Sp, z0)*pow(Sp, 1-z0), 1-z1);
    double pnew = p1/(p1+p0);
    ysample[i] = R::rbinom(1, pnew);
  }
  
  return(ysample);
}

