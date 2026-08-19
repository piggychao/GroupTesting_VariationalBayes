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


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cov_beta_cpp(
    const NumericMatrix& input_x,
    const NumericMatrix& cov_Y_vals,
    const NumericVector& Ey_1mins,
    const int psz) {
  
  int J = cov_Y_vals.nrow();
  int p = input_x.ncol();
  
  NumericMatrix cov_beta(p, p);
  
  NumericMatrix Sigma(psz, psz);
  NumericMatrix Xsub(psz, p);
  
  for (int j = 0; j < J; j++) {
    
    int offset = j * psz;
    
    //--------------------------------------------------
    // Construct Sigma_j
    //--------------------------------------------------
    
    std::fill(Sigma.begin(), Sigma.end(), 0.0);
    
    for (int i = 0; i < psz; i++)
      Sigma(i, i) = Ey_1mins[offset + i];
    
    int ct = 0;
    
    for (int ii = 1; ii < psz; ii++) {
      for (int jj = 0; jj < ii; jj++) {
        
        double val = cov_Y_vals(j, ct);
        
        Sigma(ii, jj) = val;
        Sigma(jj, ii) = val;
        
        ct++;
      }
    }
    
    //--------------------------------------------------
    // Extract X_j
    //--------------------------------------------------
    
    for (int i = 0; i < psz; i++)
      for (int k = 0; k < p; k++)
        Xsub(i, k) = input_x(offset + i, k);
    
    //--------------------------------------------------
    // cov_beta += X' Sigma X
    //--------------------------------------------------
    
    for (int a = 0; a < p; a++) {
      
      for (int b = 0; b <= a; b++) {
        
        double sum = 0.0;
        
        for (int i = 0; i < psz; i++) {
          
          for (int k = 0; k < psz; k++) {
            
            sum += Xsub(i, a) *
              Sigma(i, k) *
              Xsub(k, b);
            
          }
        }
        
        cov_beta(a, b) += sum;
        
        if (a != b)
          cov_beta(b, a) += sum;
      }
    }
  }
  
  return cov_beta;
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cal_cov_loop_cpp(
    const NumericVector& mean_beta,
    const NumericVector& mean_y,
    const NumericMatrix& input_x,
    const IntegerVector& input_z,
    const IntegerVector& input_retest_y,
    const NumericVector& input_Se,
    const NumericVector& input_Sp,
    const NumericVector& input_Se_alp_beta,
    const NumericVector& input_Sp_alp_beta,
    const IntegerMatrix& input_pool_index,
    const NumericVector& mean_z,
    const IntegerMatrix& y_mat) {
  
  int J = input_z.size();
  int size = input_pool_index.ncol();
  int p = input_x.ncol();
  int nConfig = y_mat.nrow();
  int npairs = size * (size - 1) / 2;
  
  NumericVector Ez_1mins = 1.0 - mean_z;
  
  NumericMatrix Eyy_cov(J, npairs);
  
  NumericMatrix Eyy_cov_sum11(J, npairs);
  NumericMatrix Eyy_cov_sum22(J, npairs);
  NumericMatrix Eyy_cov_sum33(J, npairs);
  
  NumericMatrix Eyy_cov_sum1(J, npairs);
  NumericMatrix Eyy_cov_sum2(J, npairs);
  NumericMatrix Eyy_cov_sum3(J, npairs);
  
  NumericMatrix cov_Se1par_sum(J, p);
  NumericMatrix cov_Sp1par_sum(J, p);
  
  NumericMatrix cov_Se2par_sum(J, p);
  NumericMatrix cov_Sp2par_sum(J, p);
  
  NumericVector Eyy_cov_sum44(J);
  NumericVector Eyy_cov_sum55(J);
  NumericVector cov_Sp2Se2_sum(J);
  
  NumericVector Eyy(npairs);
  
  //------------------------------------------------------------------
  // Main loop
  //------------------------------------------------------------------
  
  for (int j = 0; j < J; j++) {
    
    //--------------------------------------------------------------
    // Extract pool
    //--------------------------------------------------------------
    
    NumericMatrix x_sub(size, p);
    NumericVector Ey_spool(size);
    IntegerVector y_sub(size);
    NumericVector eta(size);
    
    for (int i = 0; i < size; i++) {
      
      int idx = input_pool_index(j, i) - 1;
      
      Ey_spool[i] = mean_y[idx];
      
      if (input_z[j] == 1)
        y_sub[i] = input_retest_y[idx];
      
      for (int k = 0; k < p; k++)
        x_sub(i, k) = input_x(idx, k);
    }
    
    //--------------------------------------------------------------
    // eta = X beta
    //--------------------------------------------------------------
    
    for (int i = 0; i < size; i++) {
      
      double s = 0.0;
      
      for (int k = 0; k < p; k++)
        s += x_sub(i, k) * mean_beta[k];
      
      eta[i] = s;
    }
    
    //--------------------------------------------------------------
    // Cov(Se1,beta) and Cov(Sp1,beta)
    //--------------------------------------------------------------
    
    double facSe = (input_z[j] - input_Se[0]) * Ez_1mins[j];
    double facSp = (1.0 - input_z[j] - input_Sp[0]) * (-Ez_1mins[j]);
    
    for (int k = 0; k < p; k++) {
      
      double sum = 0.0;
      
      for (int i = 0; i < size; i++)
        sum += x_sub(i, k) * Ey_spool[i];
      
      cov_Se1par_sum(j, k) = facSe * sum;
      cov_Sp1par_sum(j, k) = facSp * sum;
    }
    
    //--------------------------------------------------------------
    // Compute probabilities
    //--------------------------------------------------------------
    
    NumericVector p_sub(nConfig);
    double denom = 0.0;
    
    if (input_z[j] == 0) {
      
      for (int r = 0; r < nConfig; r++) {
        
        IntegerVector yy = y_mat(r, _);
        
        p_sub[r] = cal_prob0_cpp(
          yy,
          eta,
          input_Se_alp_beta[0],
                           input_Se_alp_beta[1],
                                            input_Sp_alp_beta[0],
                                                             input_Sp_alp_beta[1]);
        
        denom += p_sub[r];
      }
      
    } else {
      
      for (int r = 0; r < nConfig; r++) {
        
        IntegerVector yy = y_mat(r, _);
        
        p_sub[r] = cal_prob1_cpp(
          yy,
          y_sub,
          eta,
          input_Se_alp_beta[0],
                           input_Se_alp_beta[1],
                                            input_Se_alp_beta[2],
                                                             input_Se_alp_beta[3],
                                                                              input_Sp_alp_beta[0],
                                                                                               input_Sp_alp_beta[1],
                                                                                                                input_Sp_alp_beta[2],
                                                                                                                                 input_Sp_alp_beta[3]);
        
        denom += p_sub[r];
      }
    }
    
    //--------------------------------------------------------------
    // E(y_i y_j)
    //--------------------------------------------------------------
    
    int ct = 0;
    
    for (int ii = 0; ii < size - 1; ii++) {
      
      for (int jj = ii + 1; jj < size; jj++) {
        
        double num = 0.0;
        
        for (int r = 0; r < nConfig; r++)
          if (y_mat(r, ii) && y_mat(r, jj))
            num += p_sub[r];
          
          Eyy[ct] = num / denom;
          Eyy_cov(j, ct) = Eyy[ct] - Ey_spool[ii] * Ey_spool[jj];
          
          ct++;
      }
    }
    
    //--------------------------------------------------------------
    // Additional covariance terms for positive pools
    //--------------------------------------------------------------
    
    if (input_z[j] == 1) {
      
      ct = 0;
      
      for (int ii = 0; ii < size - 1; ii++) {
        
        for (int jj = ii + 1; jj < size; jj++) {
          
          double cov = Eyy_cov(j, ct);
          
          Eyy_cov_sum44[j] +=
            2.0 *
            (y_sub[ii] - input_Se[1]) *
            (y_sub[jj] - input_Se[1]) *
            cov;
          
          Eyy_cov_sum55[j] +=
            2.0 *
            (1 - y_sub[ii] - input_Sp[1]) *
            (1 - y_sub[jj] - input_Sp[1]) *
            cov;
          
          for (int k = 0; k < p; k++) {
            
            cov_Se2par_sum(j, k) +=
              ((y_sub[ii] - input_Se[1]) * x_sub(jj, k) +
              (y_sub[jj] - input_Se[1]) * x_sub(ii, k)) * cov;
            
            cov_Sp2par_sum(j, k) +=
              ((1 - y_sub[ii] - input_Sp[1]) * x_sub(jj, k) +
              (1 - y_sub[jj] - input_Sp[1]) * x_sub(ii, k)) *
              (-cov);
          }
          
          cov_Sp2Se2_sum[j] +=
            ((y_sub[ii] - input_Se[1]) *
            (1 - y_sub[jj] - input_Sp[1]) +
            (y_sub[jj] - input_Se[1]) *
            (1 - y_sub[ii] - input_Sp[1])) *
            (-cov);
          
          ct++;
        }
      }
      
      for (int ii = 0; ii < size; ii++) {
        
        cov_Sp2Se2_sum[j] +=
          ((y_sub[ii] - input_Se[1]) *
          (1 - y_sub[ii] - input_Sp[1])) *
          (-Ey_spool[ii] * (1.0 - Ey_spool[ii]));
        
        for (int k = 0; k < p; k++) {
          
          cov_Se2par_sum(j, k) +=
            (y_sub[ii] - input_Se[1]) *
            x_sub(ii, k) *
            Ey_spool[ii] *
            (1.0 - Ey_spool[ii]);
          
          cov_Sp2par_sum(j, k) +=
            (1 - y_sub[ii] - input_Sp[1]) *
            x_sub(ii, k) *
            (-Ey_spool[ii] * (1.0 - Ey_spool[ii]));
        }
      }
    }
    
    //--------------------------------------------------------------
    // x-sub products
    //--------------------------------------------------------------
    
    ct = 0;
    
    for (int ii = 0; ii < size - 1; ii++) {
      
      for (int jj = ii + 1; jj < size; jj++) {
        
        double cov = Eyy_cov(j, ct);
        
        Eyy_cov_sum11(j, ct) =
          2.0 * x_sub(ii, 0) * x_sub(jj, 0) * cov;
        
        Eyy_cov_sum22(j, ct) =
          2.0 * x_sub(ii, 1) * x_sub(jj, 1) * cov;
        
        Eyy_cov_sum33(j, ct) =
          2.0 * x_sub(ii, 2) * x_sub(jj, 2) * cov;
        
        Eyy_cov_sum1(j, ct) =
          (x_sub(ii, 1) + x_sub(jj, 1)) * cov;
        
        Eyy_cov_sum2(j, ct) =
          (x_sub(ii, 2) + x_sub(jj, 2)) * cov;
        
        Eyy_cov_sum3(j, ct) =
          (x_sub(ii, 1) * x_sub(jj, 2) +
          x_sub(jj, 1) * x_sub(ii, 2)) * cov;
        
        ct++;
      }
    }
  }
  
  return List::create(
    _["Eyy_cov"] = Eyy_cov,
    _["Eyy_cov_sum11"] = Eyy_cov_sum11,
    _["Eyy_cov_sum22"] = Eyy_cov_sum22,
    _["Eyy_cov_sum33"] = Eyy_cov_sum33,
    _["Eyy_cov_sum44"] = Eyy_cov_sum44,
    _["Eyy_cov_sum55"] = Eyy_cov_sum55,
    _["Eyy_cov_sum1"] = Eyy_cov_sum1,
    _["Eyy_cov_sum2"] = Eyy_cov_sum2,
    _["Eyy_cov_sum3"] = Eyy_cov_sum3,
    _["cov_Se1par_sum"] = cov_Se1par_sum,
    _["cov_Sp1par_sum"] = cov_Sp1par_sum,
    _["cov_Se2par_sum"] = cov_Se2par_sum,
    _["cov_Sp2par_sum"] = cov_Sp2par_sum,
    _["cov_Sp2Se2_sum"] = cov_Sp2Se2_sum
  );
}