#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cal_prob0_cpp(const IntegerVector& yvals,
                     const NumericVector& x_abc_input_beta,
                     double Se1_alpha,
                     double Se1_beta,
                     double Sp1_alpha,
                     double Sp1_beta)
{
  int pool_abc = 0;
  double eta = 0.0;
  
  int n = yvals.size();
  
  for (int i = 0; i < n; i++) {
    pool_abc += yvals[i];
    eta += yvals[i] * x_abc_input_beta[i];
  }
  
  double log_1minSe1 =
    R::digamma(Se1_beta) -
    R::digamma(Se1_alpha + Se1_beta);
  
  double log_Sp1 =
    R::digamma(Sp1_alpha) -
    R::digamma(Sp1_alpha + Sp1_beta);
  
  double logp;
  
  if (pool_abc > 0)
    logp = log_1minSe1 + eta;
  else
    logp = log_Sp1 + eta;
  
  return std::exp(logp);
}
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector VB_ytilde0_cpp(const NumericMatrix& x_sub,
                             const NumericVector& mu_beta,
                             const IntegerMatrix& y_mat,
                             double Se_alpha,
                             double Se_beta,
                             double Sp_alpha,
                             double Sp_beta) {
  
  int pool_size = x_sub.nrow();
  int nConfig   = y_mat.nrow();
  int p         = x_sub.ncol();
  
  //-------------------------------------------------------
  // Compute x_sub %*% mu_beta
  //-------------------------------------------------------
  
  NumericVector eta(pool_size);
  
  for (int i = 0; i < pool_size; i++) {
    double s = 0.0;
    for (int j = 0; j < p; j++)
      s += x_sub(i,j) * mu_beta[j];
    eta[i] = s;
  }
  
  //-------------------------------------------------------
  // Compute probabilities
  //-------------------------------------------------------
  
  NumericVector p_sub(nConfig);
  
  double denom = 0.0;
  
  for (int i = 0; i < nConfig; i++) {
    
    IntegerVector y = y_mat(i, _);
    
    p_sub[i] = cal_prob0_cpp(y,
                             eta,
                             Se_alpha,
                             Se_beta,
                             Sp_alpha,
                             Sp_beta);
    
    denom += p_sub[i];
  }
  
  //-------------------------------------------------------
  // Posterior probabilities
  //-------------------------------------------------------
  
  NumericVector out(pool_size + 1);
  
  for (int ind = 0; ind < pool_size; ind++) {
    
    double num = 0.0;
    
    for (int config = 0; config < nConfig; config++) {
      if (y_mat(config, ind) == 1)
        num += p_sub[config];
    }
    
    out[ind] = num / denom;
  }
  
  out[pool_size] = p_sub[0] / denom;
  
  return out;
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cal_prob1_cpp(const IntegerVector& yvals,
                     const IntegerVector& y_abc,
                     const NumericVector& x_abc_input_beta,
                     double Se1_alpha,
                     double Se1_beta,
                     double Se2_alpha,
                     double Se2_beta,
                     double Sp1_alpha,
                     double Sp1_beta,
                     double Sp2_alpha,
                     double Sp2_beta)
{
  int n = yvals.size();
  
  bool positive_pool = false;
  
  double eta = 0.0;
  double test2 = 0.0;
  
  for (int i = 0; i < n; i++) {
    
    if (yvals[i])
      positive_pool = true;
    
    eta += yvals[i] * x_abc_input_beta[i];
    
    if (yvals[i]) {
      
      if (y_abc[i])
        test2 += R::digamma(Se2_alpha) -
          R::digamma(Se2_alpha + Se2_beta);
      else
        test2 += R::digamma(Se2_beta) -
          R::digamma(Se2_alpha + Se2_beta);
      
    } else {
      
      if (y_abc[i])
        test2 += R::digamma(Sp2_beta) -
          R::digamma(Sp2_alpha + Sp2_beta);
      else
        test2 += R::digamma(Sp2_alpha) -
          R::digamma(Sp2_alpha + Sp2_beta);
    }
  }
  
  double log1;
  
  if (positive_pool)
    log1 = R::digamma(Se1_alpha) -
      R::digamma(Se1_alpha + Se1_beta);
  else
    log1 = R::digamma(Sp1_beta) -
      R::digamma(Sp1_alpha + Sp1_beta);
  
  return std::exp(log1 + test2 + eta);
}

// [[Rcpp::export]]
Rcpp::NumericVector VB_ytilde1_cpp(
    const Rcpp::IntegerVector& input_z_sub,
    const Rcpp::NumericMatrix& input_x_sub,
    const Rcpp::NumericVector& input_mu_beta,
    const Rcpp::IntegerMatrix& y_mat,
    double Se1_alpha,
    double Se1_beta,
    double Se2_alpha,
    double Se2_beta,
    double Sp1_alpha,
    double Sp1_beta,
    double Sp2_alpha,
    double Sp2_beta)
{
  using namespace Rcpp;
  
  int pool_size = input_x_sub.nrow();
  int p = input_x_sub.ncol();
  int nConfig = y_mat.nrow();
  
  //-------------------------------------------------------
  // Compute x_sub %*% beta
  //-------------------------------------------------------
  
  NumericVector eta(pool_size);
  
  for (int i = 0; i < pool_size; i++) {
    double s = 0.0;
    for (int j = 0; j < p; j++)
      s += input_x_sub(i, j) * input_mu_beta[j];
    
    eta[i] = s;
  }
  
  //-------------------------------------------------------
  // Compute probabilities
  //-------------------------------------------------------
  
  NumericVector p_sub(nConfig);
  
  double denom = 0.0;
  
  for (int config = 0; config < nConfig; config++) {
    
    IntegerVector y = y_mat(config, _);
    
    p_sub[config] =
      cal_prob1_cpp(
        y,
        input_z_sub,
        eta,
        Se1_alpha,
        Se1_beta,
        Se2_alpha,
        Se2_beta,
        Sp1_alpha,
        Sp1_beta,
        Sp2_alpha,
        Sp2_beta);
    
    denom += p_sub[config];
  }
  
  //-------------------------------------------------------
  // Posterior probabilities
  //-------------------------------------------------------
  
  NumericVector out(pool_size + 1);
  
  for (int ind = 0; ind < pool_size; ind++) {
    
    double num = 0.0;
    
    for (int config = 0; config < nConfig; config++) {
      
      if (y_mat(config, ind) == 1)
        num += p_sub[config];
    }
    
    out[ind] = num / denom;
  }
  
  out[pool_size] = p_sub[0] / denom;
  
  return out;
}

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List cov_sums_loop_cpp(
    NumericMatrix input_x,
    IntegerVector input_z,
    IntegerVector input_retest_y,
    NumericVector mean_y,
    NumericVector mean_beta,
    NumericVector input_Se,
    NumericVector input_Sp,
    NumericVector input_Se_alp_beta,
    NumericVector input_Sp_alp_beta,
    IntegerMatrix input_pool_index,
    IntegerMatrix y_mat,
    NumericVector Ez_1mins
){
  
  int J = input_pool_index.nrow();
  int size = input_pool_index.ncol();
  int p = input_x.ncol();
  
  int nPairs = size*(size-1)/2;
  
  
  NumericMatrix cov_Se1par_sum(J,p);
  NumericMatrix cov_Sp1par_sum(J,p);
  
  NumericMatrix Eyy_cov(J,nPairs);
  
  NumericVector Eyy(nPairs);
  
  NumericVector Eyy_cov_sum44(J);
  NumericVector Eyy_cov_sum55(J);
  
  NumericVector cov_Se2par_sum1(J);
  NumericVector cov_Se2par_sum2(J);
  NumericVector cov_Se2par_sum3(J);
  
  NumericVector cov_Sp2par_sum1(J);
  NumericVector cov_Sp2par_sum2(J);
  NumericVector cov_Sp2par_sum3(J);
  
  NumericVector cov_Sp2Se2_sum(J);
  
  
  NumericMatrix Eyy_cov_sum11(J,nPairs);
  NumericMatrix Eyy_cov_sum22(J,nPairs);
  NumericMatrix Eyy_cov_sum33(J,nPairs);
  
  NumericMatrix Eyy_cov_sum1(J,nPairs);
  NumericMatrix Eyy_cov_sum2(J,nPairs);
  NumericMatrix Eyy_cov_sum3(J,nPairs);
  
  
  
  for(int j=0; j<J; j++){
    
    // -----------------------------
    // Extract pool information
    // -----------------------------
    
    IntegerVector pool_index = input_pool_index(j,_);
    
    
    NumericMatrix x_sub(size,p);
    
    NumericVector Ey_spool(size);
    
    
    for(int i=0;i<size;i++){
      
      int ind = pool_index[i]-1; // R -> C index
      
      Ey_spool[i] = mean_y[ind];
      
      for(int k=0;k<p;k++)
        x_sub(i,k)=input_x(ind,k);
      
    }
    
    
    // x_sub %*% beta
    
    NumericVector x_sub_mean_beta(size);
    
    for(int i=0;i<size;i++){
      
      double tmp=0;
      
      for(int k=0;k<p;k++)
        tmp += x_sub(i,k)*mean_beta[k];
      
      x_sub_mean_beta[i]=tmp;
    }
    
    
    
    // -----------------------------
    // Cov Se1 / Sp1 with beta
    // -----------------------------
    
    for(int k=0;k<p;k++){
      
      double s1=0;
      double s2=0;
      
      for(int i=0;i<size;i++){
        
        s1 += x_sub(i,k) *
          Ey_spool[i] *
          Ez_1mins[j];
        
        s2 += x_sub(i,k) *
          Ey_spool[i] *
          (-Ez_1mins[j]);
        
      }
      
      
      cov_Se1par_sum(j,k)=
        (input_z[j]-input_Se[0])*s1;
      
      
      cov_Sp1par_sum(j,k)=
        (1-input_z[j]-input_Sp[0])*s2;
      
    }
    
    
    
    // ======================================================
    // Calculate Eyy covariance
    // ======================================================
    
    
    if(input_z[j]==0){
      
      
      NumericVector p_sub(y_mat.nrow());
      
      double denom=0;
      
      
      for(int config=0; config<y_mat.nrow(); config++){
        
        IntegerVector row=y_mat(config,_);
        
        
        p_sub[config]=cal_prob0_cpp(
          row,
          x_sub_mean_beta,
          input_Se_alp_beta[0],
                           input_Se_alp_beta[1],
                                            input_Sp_alp_beta[0],
                                                             input_Sp_alp_beta[1]
        );
        
        denom+=p_sub[config];
        
      }
      
      
      
      int ct=0;
      
      
      for(int ii=0;ii<size-1;ii++){
        
        for(int jj=ii+1;jj<size;jj++){
          
          
          double numerator=0;
          
          
          for(int config=0;config<y_mat.nrow();config++){
            
            if(y_mat(config,ii)==1 &&
               y_mat(config,jj)==1)
              
              numerator += p_sub[config];
            
          }
          
          
          Eyy[ct]=numerator/denom;
          
          ct++;
          
        }
      }
      
      
      
      ct=0;
      
      for(int ii=0;ii<size-1;ii++){
        
        for(int jj=ii+1;jj<size;jj++){
          
          Eyy_cov(j,ct)=
            Eyy[ct] -
            Ey_spool[ii]*Ey_spool[jj];
          
          ct++;
          
        }
      }
      
    }
    
    
    
    
    // ======================================================
    // Positive pools
    // ======================================================
    
    if(input_z[j]==1){
      
      
      IntegerVector y_sub(size);
      
      
      for(int i=0;i<size;i++)
        y_sub[i]=input_retest_y[pool_index[i]-1];
      
      
      
      NumericVector p_sub(y_mat.nrow());
      
      double denom=0;
      
      
      for(int config=0;config<y_mat.nrow();config++){
        
        IntegerVector row=y_mat(config,_);
        
        
        p_sub[config]=cal_prob1_cpp(
          row,
          y_sub,
          x_sub_mean_beta,
          input_Se_alp_beta[0],
                           input_Se_alp_beta[1],
                                            input_Se_alp_beta[2],
                                                             input_Se_alp_beta[3],
                                                                              input_Sp_alp_beta[0],
                                                                                               input_Sp_alp_beta[1],
                                                                                                                input_Sp_alp_beta[2],
                                                                                                                                 input_Sp_alp_beta[3]
        );
        
        
        denom+=p_sub[config];
        
      }
      
      
      
      int ct=0;
      
      
      for(int ii=0;ii<size-1;ii++){
        
        for(int jj=ii+1;jj<size;jj++){
          
          
          double numerator=0;
          
          
          for(int config=0;config<y_mat.nrow();config++){
            
            if(y_mat(config,ii)==1 &&
               y_mat(config,jj)==1)
              
              numerator+=p_sub[config];
            
          }
          
          
          Eyy[ct]=numerator/denom;
          
          ct++;
          
        }
      }
      
      
      
      ct=0;
      
      
      for(int ii=0;ii<size-1;ii++){
        
        for(int jj=ii+1;jj<size;jj++){
          
          
          double cov=Eyy[ct]-
            Ey_spool[ii]*Ey_spool[jj];
          
          
          Eyy_cov(j,ct)=cov;
          
          
          
          double sei =
            y_sub[ii]-input_Se[1];
          
          double sej =
            y_sub[jj]-input_Se[1];
          
          
          Eyy_cov_sum44[j]+=2*sei*sej*cov;
          
          
          double spi =
            1-y_sub[ii]-input_Sp[1];
          
          double spj =
            1-y_sub[jj]-input_Sp[1];
          
          
          Eyy_cov_sum55[j]+=2*spi*spj*cov;
          
          
          cov_Se2par_sum1[j]+=(sei+sej)*cov;
          
          
          cov_Se2par_sum2[j]+=
            (sei*x_sub(jj,1)+
            sej*x_sub(ii,1))*cov;
          
          
          cov_Se2par_sum3[j]+=
            (sei*x_sub(jj,2)+
            sej*x_sub(ii,2))*cov;
          
          
          
          cov_Sp2par_sum1[j]+=
            (spi+spj)*(-cov);
          
          
          cov_Sp2par_sum2[j]+=
            (spi*x_sub(jj,1)+
            spj*x_sub(ii,1))*(-cov);
          
          
          cov_Sp2par_sum3[j]+=
            (spi*x_sub(jj,2)+
            spj*x_sub(ii,2))*(-cov);
          
          
          
          cov_Sp2Se2_sum[j]+=
            ((Ey_spool[ii]-input_Se[1])*
            (1-Ey_spool[jj]-input_Sp[1])
               +
                 (Ey_spool[jj]-input_Se[1])*
                 (1-Ey_spool[ii]-input_Sp[1]))
            *
              (-cov);
          
          
          ct++;
          
        }
      }
      
      
      
      for(int i=0;i<size;i++){
        
        
        double tmp1 =
          (y_sub[i]-input_Se[1])*
          Ey_spool[i]*
          (1-Ey_spool[i]);
        
        
        cov_Se2par_sum1[j]+=tmp1;
        
        
        cov_Se2par_sum2[j]+=
          (y_sub[i]-input_Se[1])*
          x_sub(i,1)*
          Ey_spool[i]*
          (1-Ey_spool[i]);
        
        
        cov_Se2par_sum3[j]+=
          (y_sub[i]-input_Se[1])*
          x_sub(i,2)*
          Ey_spool[i]*
          (1-Ey_spool[i]);
        
      }
      
    }
    
    
    
    // ======================================================
    // Beta-beta covariance terms
    // ======================================================
    
    
    int ct=0;
    
    for(int ii=0;ii<size-1;ii++){
      
      for(int jj=ii+1;jj<size;jj++){
        
        
        double cov=Eyy_cov(j,ct);
        
        
        Eyy_cov_sum11(j,ct)=
          (x_sub(ii,0)*x_sub(jj,0)*2)*cov;
        
        
        Eyy_cov_sum22(j,ct)=
          (x_sub(ii,1)*x_sub(jj,1)*2)*cov;
        
        
        Eyy_cov_sum33(j,ct)=
          (x_sub(ii,2)*x_sub(jj,2)*2)*cov;
        
        
        
        Eyy_cov_sum1(j,ct)=
          (x_sub(ii,1)+x_sub(jj,1))*cov;
        
        
        Eyy_cov_sum2(j,ct)=
          (x_sub(ii,2)+x_sub(jj,2))*cov;
        
        
        Eyy_cov_sum3(j,ct)=
          (x_sub(ii,1)*x_sub(jj,2)+
          x_sub(jj,1)*x_sub(ii,2))*cov;
        
        
        ct++;
        
      }
    }
    
    
  }
  
  
  return List::create(
    _["cov_Se1par_sum"]=cov_Se1par_sum,
    _["cov_Sp1par_sum"]=cov_Sp1par_sum,
    _["Eyy_cov"]=Eyy_cov,
    _["Eyy_cov_sum11"]=Eyy_cov_sum11,
    _["Eyy_cov_sum22"]=Eyy_cov_sum22,
    _["Eyy_cov_sum33"]=Eyy_cov_sum33,
    _["Eyy_cov_sum1"]=Eyy_cov_sum1,
    _["Eyy_cov_sum2"]=Eyy_cov_sum2,
    _["Eyy_cov_sum3"]=Eyy_cov_sum3,
    _["Eyy_cov_sum44"]=Eyy_cov_sum44,
    _["Eyy_cov_sum55"]=Eyy_cov_sum55,
    _["cov_Se2par_sum1"]=cov_Se2par_sum1,
    _["cov_Se2par_sum2"]=cov_Se2par_sum2,
    _["cov_Se2par_sum3"]=cov_Se2par_sum3,
    _["cov_Sp2par_sum1"]=cov_Sp2par_sum1,
    _["cov_Sp2par_sum2"]=cov_Sp2par_sum2,
    _["cov_Sp2par_sum3"]=cov_Sp2par_sum3,
    _["cov_Sp2Se2_sum"]=cov_Sp2Se2_sum
  );
  
}

#include <Rcpp.h>

using namespace Rcpp;

inline double cal_prob0_cpp(const NumericMatrix& y_mat,
                            int row_idx,
                            const std::vector<double>& xbeta,
                            double Se1_alpha, double Se1_beta,
                            double Sp1_alpha, double Sp1_beta) {
  const int ncol = y_mat.ncol();
  
  double pool_abc = 0.0;
  double x_term = 0.0;
  
  for (int i = 0; i < ncol; ++i) {
    const double y = y_mat(row_idx, i);
    pool_abc += y;
    x_term += y * xbeta[i];
  }
  
  const double log_1minSe1 = R::digamma(Se1_beta) - R::digamma(Se1_alpha + Se1_beta);
  const double log_Sp1     = R::digamma(Sp1_alpha) - R::digamma(Sp1_alpha + Sp1_beta);
  
  const double ll = (pool_abc > 0.0 ? log_1minSe1 : log_Sp1) + x_term;
  return std::exp(ll);
}

inline double cal_prob1_cpp(const NumericMatrix& y_mat,
                            int row_idx,
                            const std::vector<double>& y_abc,
                            const std::vector<double>& xbeta,
                            double Se1_alpha, double Se1_beta,
                            double Se2_alpha, double Se2_beta,
                            double Sp1_alpha, double Sp1_beta,
                            double Sp2_alpha, double Sp2_beta) {
  const int ncol = y_mat.ncol();
  
  double pool_abc = 0.0;
  double term_y = 0.0;
  double term_n = 0.0;
  double x_term = 0.0;
  
  const double log_Se1     = R::digamma(Se1_alpha) - R::digamma(Se1_alpha + Se1_beta);
  const double log_1minSp1 = R::digamma(Sp1_beta)  - R::digamma(Sp1_alpha + Sp1_beta);
  const double log_Se2     = R::digamma(Se2_alpha) - R::digamma(Se2_alpha + Se2_beta);
  const double log_1minSe2 = R::digamma(Se2_beta)  - R::digamma(Se2_alpha + Se2_beta);
  const double log_Sp2     = R::digamma(Sp2_alpha) - R::digamma(Sp2_alpha + Sp2_beta);
  const double log_1minSp2 = R::digamma(Sp2_beta)  - R::digamma(Sp2_alpha + Sp2_beta);
  
  for (int i = 0; i < ncol; ++i) {
    const double y = y_mat(row_idx, i);
    const double ya = y_abc[i];
    
    pool_abc += y;
    term_y += y * (ya * log_Se2 + (1.0 - ya) * log_1minSe2);
    term_n += (1.0 - y) * ((1.0 - ya) * log_Sp2 + ya * log_1minSp2);
    x_term += y * xbeta[i];
  }
  
  const double ll = (pool_abc > 0.0 ? log_Se1 : log_1minSp1) + term_y + term_n + x_term;
  return std::exp(ll);
}
#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

inline double log_prob0_cpp(const NumericVector& yrow,
                            const NumericVector& xbeta,
                            double log_1minSe1,
                            double log_Sp1) {
  const int n = yrow.size();
  double pool_abc = 0.0;
  double linpred  = 0.0;
  
  for (int k = 0; k < n; ++k) {
    const double y = yrow[k];
    pool_abc += y;
    linpred  += y * xbeta[k];
  }
  
  return ((pool_abc > 0.0) ? log_1minSe1 : log_Sp1) + linpred;
}

inline double log_prob1_cpp(const NumericVector& yrow,
                            const NumericVector& y_sub,
                            const NumericVector& xbeta,
                            double log_Se1,
                            double log_1minSp1,
                            double log_Se2,
                            double log_1minSe2,
                            double log_Sp2,
                            double log_1minSp2) {
  const int n = yrow.size();
  double pool_abc = 0.0;
  double linpred  = 0.0;
  double term2    = 0.0;
  double term3    = 0.0;
  
  for (int k = 0; k < n; ++k) {
    const double y  = yrow[k];
    const double ys = y_sub[k];
    
    pool_abc += y;
    linpred  += y * xbeta[k];
    
    term2 += y * (ys * log_Se2 + (1.0 - ys) * log_1minSe2);
    term3 += (1.0 - y) * ((1.0 - ys) * log_Sp2 + ys * log_1minSp2);
  }
  
  return ((pool_abc > 0.0) ? log_Se1 : log_1minSp1) + term2 + term3 + linpred;
}

// [[Rcpp::export]]
NumericVector cal_cov_sums_unknw_cpp(const NumericVector& mean_beta,
                                     const NumericVector& mean_y,
                                     const NumericMatrix& input_x,
                                     const NumericVector& input_z,
                                     const NumericVector& input_retest_y,
                                     const NumericVector& input_Se,
                                     const NumericVector& input_Sp,
                                     const NumericVector& input_Se_alp_beta,
                                     const NumericVector& input_Sp_alp_beta,
                                     const NumericMatrix& input_pool_index,
                                     const NumericVector& mean_z,
                                     const NumericMatrix& y_mat) {
  
  const int J    = input_z.size();
  const int size = input_pool_index.ncol();
  const int p    = input_x.ncol();
  const int nconf = y_mat.nrow();
  
  if (p < 3) {
    stop("input_x must have at least 3 columns (intercept + x1 + x2), as in the original R code.");
  }
  if (input_pool_index.nrow() != J) {
    stop("input_pool_index must have J rows, where J = length(input_z).");
  }
  if (mean_z.size() != J) {
    stop("mean_z must have the same length as input_z.");
  }
  
  // Precompute Beta-log terms once (same values for every pool).
  const double log_1minSe1 = R::digamma(input_Se_alp_beta[1]) -
    R::digamma(input_Se_alp_beta[0] + input_Se_alp_beta[1]);
  const double log_Se1     = R::digamma(input_Se_alp_beta[0]) -
    R::digamma(input_Se_alp_beta[0] + input_Se_alp_beta[1]);
  const double log_1minSp1 = R::digamma(input_Sp_alp_beta[1]) -
    R::digamma(input_Sp_alp_beta[0] + input_Sp_alp_beta[1]);
  const double log_Sp1     = R::digamma(input_Sp_alp_beta[0]) -
    R::digamma(input_Sp_alp_beta[0] + input_Sp_alp_beta[1]);
  
  const double log_Se2     = R::digamma(input_Se_alp_beta[2]) -
    R::digamma(input_Se_alp_beta[2] + input_Se_alp_beta[3]);
  const double log_1minSe2 = R::digamma(input_Se_alp_beta[3]) -
    R::digamma(input_Se_alp_beta[2] + input_Se_alp_beta[3]);
  const double log_Sp2     = R::digamma(input_Sp_alp_beta[2]) -
    R::digamma(input_Sp_alp_beta[2] + input_Sp_alp_beta[3]);
  const double log_1minSp2 = R::digamma(input_Sp_alp_beta[3]) -
    R::digamma(input_Sp_alp_beta[2] + input_Sp_alp_beta[3]);
  
  const int pair_count = size * (size - 1) / 2;
  
  double sum_Eyy_cov_sum11 = 0.0;
  double sum_Eyy_cov_sum22 = 0.0;
  double sum_Eyy_cov_sum33 = 0.0;
  
  double sum_Eyy_cov_sum44 = 0.0;
  double sum_Eyy_cov_sum55 = 0.0;
  
  double sum_Eyy_cov_sum1  = 0.0;
  double sum_Eyy_cov_sum2  = 0.0;
  double sum_Eyy_cov_sum3  = 0.0;
  
  std::vector<double> cov_Se1par_sum(p, 0.0);
  std::vector<double> cov_Sp1par_sum(p, 0.0);
  
  std::vector<double> cov_Se2par_sum1(J, 0.0);
  std::vector<double> cov_Se2par_sum2(J, 0.0);
  std::vector<double> cov_Se2par_sum3(J, 0.0);
  
  std::vector<double> cov_Sp2par_sum1(J, 0.0);
  std::vector<double> cov_Sp2par_sum2(J, 0.0);
  std::vector<double> cov_Sp2par_sum3(J, 0.0);
  
  std::vector<double> cov_Sp2Se2_sum(J, 0.0);
  
  for (int j = 0; j < J; ++j) {
    const double z_j   = input_z[j];
    const double Se1   = input_Se[0];
    const double Se2   = input_Se[1];
    const double Sp1   = input_Sp[0];
    const double Sp2   = input_Sp[1];
    const double Ez1m  = 1.0 - mean_z[j];
    const double Ez    = mean_z[j];
    
    NumericVector Ey_spool(size);
    NumericVector y_sub(size);
    NumericVector xbeta(size);
    
    for (int ii = 0; ii < size; ++ii) {
      const int idx = static_cast<int>(input_pool_index(j, ii)) - 1;
      Ey_spool[ii] = mean_y[idx];
      
      double xb = 0.0;
      for (int k = 0; k < p; ++k) {
        xb += input_x(idx, k) * mean_beta[k];
      }
      xbeta[ii] = xb;
      
      if (z_j == 1.0) {
        y_sub[ii] = input_retest_y[idx];
      }
    }
    
    // cov_Se1par_sum and cov_Sp1par_sum
    for (int k = 0; k < p; ++k) {
      double c1 = 0.0;
      double c2 = 0.0;
      for (int ii = 0; ii < size; ++ii) {
        const int idx = static_cast<int>(input_pool_index(j, ii)) - 1;
        c1 += input_x(idx, k) * Ey_spool[ii];
        c2 += input_x(idx, k) * Ey_spool[ii];
      }
      cov_Se1par_sum[k] += (z_j - Se1) * Ez1m * c1;
      cov_Sp1par_sum[k] += (1.0 - z_j - Sp1) * (-Ez1m) * c2;
    }
    
    // Compute normalized p_sub over all y.mat configurations using log-sum-exp
    std::vector<double> logw(nconf);
    for (int r = 0; r < nconf; ++r) {
      NumericVector yrow = y_mat(r, _);
      
      if (z_j == 0.0) {
        logw[r] = log_prob0_cpp(yrow, xbeta, log_1minSe1, log_Sp1);
      } else {
        logw[r] = log_prob1_cpp(yrow, y_sub, xbeta,
                                log_Se1, log_1minSp1, log_Se2, log_1minSe2,
                                log_Sp2, log_1minSp2);
      }
    }
    
    const double max_logw = *std::max_element(logw.begin(), logw.end());
    std::vector<double> w(nconf);
    double denom = 0.0;
    for (int r = 0; r < nconf; ++r) {
      w[r] = std::exp(logw[r] - max_logw);
      denom += w[r];
    }
    if (denom == 0.0 || !R_finite(denom)) {
      stop("Numerical underflow/overflow in cal_cov_sums_unknw_cpp: denominator is invalid.");
    }
    
    // Pairwise expectations and covariance sums
    int ct = 0;
    for (int ii = 0; ii < size - 1; ++ii) {
      for (int jj = ii + 1; jj < size; ++jj) {
        double num_11 = 0.0;
        for (int r = 0; r < nconf; ++r) {
          if (y_mat(r, ii) == 1.0 && y_mat(r, jj) == 1.0) {
            num_11 += w[r];
          }
        }
        
        const double Eyy = num_11 / denom;
        const double cov = Eyy - Ey_spool[ii] * Ey_spool[jj];
        
        // These reproduce x_sub_multi[1,], [2,], [3,] from the original R code
        sum_Eyy_cov_sum11 += 2.0 * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 0) *
          input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 0) * cov;
        sum_Eyy_cov_sum22 += 2.0 * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 1) *
          input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 1) * cov;
        sum_Eyy_cov_sum33 += 2.0 * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 2) *
          input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 2) * cov;
        
        sum_Eyy_cov_sum1 += (input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 1) +
          input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 1)) * cov;
        sum_Eyy_cov_sum2 += (input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 2) +
          input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 2)) * cov;
        sum_Eyy_cov_sum3 += (input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 1) *
          input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 2) +
          input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 1) *
          input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 2)) * cov;
        
        if (z_j == 1.0) {
          const double y_i = y_sub[ii];
          const double y_j = y_sub[jj];
          
          sum_Eyy_cov_sum44 += 2.0 * (y_i - Se2) * (y_j - Se2) * cov;
          sum_Eyy_cov_sum55 += 2.0 * (1.0 - y_i - Sp2) * (1.0 - y_j - Sp2) * cov;
          
          cov_Se2par_sum1[j] += ((y_i - Se2) + (y_j - Se2)) * cov;
          cov_Se2par_sum2[j] += ((y_i - Se2) * input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 1) +
            (y_j - Se2) * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 1)) * cov;
          cov_Se2par_sum3[j] += ((y_i - Se2) * input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 2) +
            (y_j - Se2) * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 2)) * cov;
          
          cov_Sp2par_sum1[j] += ((1.0 - y_i - Sp2) + (1.0 - y_j - Sp2)) * (-cov);
          cov_Sp2par_sum2[j] += ((1.0 - y_i - Sp2) * input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 1) +
            (1.0 - y_j - Sp2) * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 1)) * (-cov);
          cov_Sp2par_sum3[j] += ((1.0 - y_i - Sp2) * input_x(static_cast<int>(input_pool_index(j, jj)) - 1, 2) +
            (1.0 - y_j - Sp2) * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 2)) * (-cov);
          
          cov_Sp2Se2_sum[j] += ((y_i - Se2) * (1.0 - y_j - Sp2) +
            (y_j - Se2) * (1.0 - y_i - Sp2)) * (-cov);
        }
        
        ++ct;
      }
    }
    
    if (z_j == 1.0) {
      for (int ii = 0; ii < size; ++ii) {
        const double v = Ey_spool[ii] * (1.0 - Ey_spool[ii]);
        cov_Sp2Se2_sum[j] += ((y_sub[ii] - Se2) * (1.0 - y_sub[ii] - Sp2)) * (-v);
        
        cov_Se2par_sum1[j] += (y_sub[ii] - Se2) * v;
        cov_Se2par_sum2[j] += (y_sub[ii] - Se2) * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 1) * v;
        cov_Se2par_sum3[j] += (y_sub[ii] - Se2) * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 2) * v;
        
        cov_Sp2par_sum1[j] += (1.0 - y_sub[ii] - Sp2) * (-v);
        cov_Sp2par_sum2[j] += (1.0 - y_sub[ii] - Sp2) * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 1) * (-v);
        cov_Sp2par_sum3[j] += (1.0 - y_sub[ii] - Sp2) * input_x(static_cast<int>(input_pool_index(j, ii)) - 1, 2) * (-v);
      }
    }
  }
  
  const int out_len = 2 * p + 15;
  NumericVector out(out_len);
  int pos = 0;
  
  out[pos++] = sum_Eyy_cov_sum11;
  out[pos++] = sum_Eyy_cov_sum22;
  out[pos++] = sum_Eyy_cov_sum33;
  
  out[pos++] = sum_Eyy_cov_sum44;
  out[pos++] = sum_Eyy_cov_sum55;
  
  out[pos++] = sum_Eyy_cov_sum1;
  out[pos++] = sum_Eyy_cov_sum2;
  out[pos++] = sum_Eyy_cov_sum3;
  
  for (int k = 0; k < p; ++k) out[pos++] = cov_Se1par_sum[k];
  for (int k = 0; k < p; ++k) out[pos++] = cov_Sp1par_sum[k];
  
  out[pos++] = std::accumulate(cov_Se2par_sum1.begin(), cov_Se2par_sum1.end(), 0.0);
  out[pos++] = std::accumulate(cov_Se2par_sum2.begin(), cov_Se2par_sum2.end(), 0.0);
  out[pos++] = std::accumulate(cov_Se2par_sum3.begin(), cov_Se2par_sum3.end(), 0.0);
  
  out[pos++] = std::accumulate(cov_Sp2par_sum1.begin(), cov_Sp2par_sum1.end(), 0.0);
  out[pos++] = std::accumulate(cov_Sp2par_sum2.begin(), cov_Sp2par_sum2.end(), 0.0);
  out[pos++] = std::accumulate(cov_Sp2par_sum3.begin(), cov_Sp2par_sum3.end(), 0.0);
  
  out[pos++] = std::accumulate(cov_Sp2Se2_sum.begin(), cov_Sp2Se2_sum.end(), 0.0);
  
  return out;
}
