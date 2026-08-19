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

// pooi = pool index
// pind: current individual probability 
// z = pool response

// [[Rcpp::export]]
NumericVector sampletildey1(NumericVector pooli, NumericVector pind, 
                     NumericVector z, NumericVector yinit, 
                     double Se, double Sp, int size){
  int nsample = yinit.length();
  
  for(int i=0; i<nsample; ++i){
    int id = pooli[i]; // get pool index for individual i
    int startid = id*size;
    int endid = id*size+(size-1);
    IntegerVector poolind = seq(startid, endid); // use pool-individual ids to find subset pool
    double pind0 = pind[i]; // probability of individual i
    double z0 = z[id]; // get the corresponding pool response
    double p1 = pind0 * pow(Se, z0) * pow(1-Se, (1-z0));
    NumericVector spool = yinit[poolind]; // find corresponding subset pool
    int yind0 = i - size*id; // get individual index in subset pool
    spool.erase(spool.begin() + yind0); // removing y_i from subset pool
    double spools = std::accumulate(spool.begin(), spool.end(), 0);
    double z1;
    if (spools > 0) {z1=1;} else {z1=0;} // pool response for new subset pool
    double p0 = (1-pind0) * pow( pow(Se, z0)*pow(1-Se, 1-z0), z1) *
      pow( pow(1-Sp, z0)*pow(Sp, 1-z0), 1-z1);
    double pnew = p1/(p1+p0);
    yinit[i] = R::rbinom(1, pnew);
  }
  
  return(yinit);
}


