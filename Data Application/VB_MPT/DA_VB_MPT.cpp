#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector VB_ytilde0_cpp(const NumericMatrix& x_sub,
                             const NumericVector& mu_beta,
                             const NumericMatrix& y_mat,
                             const double Se,
                             const double Sp) {
  
  int pool_size = x_sub.nrow();
  int n_config = y_mat.nrow();
  
  // x_sub %*% mu_beta
  NumericVector eta(pool_size);
  for(int i = 0; i < pool_size; i++) {
    double tmp = 0.0;
    for(int j = 0; j < x_sub.ncol(); j++)
      tmp += x_sub(i,j) * mu_beta[j];
      eta[i] = tmp;
  }
  
  NumericVector p_sub(n_config);
  
  //-------------------------------------------------------
    // Calculate p_sub
  //-------------------------------------------------------
    
    for(int config = 0; config < n_config; config++) {
      
      int pool_sum = 0;
      double linear = 0.0;
      
      for(int j = 0; j < pool_size; j++) {
        double y = y_mat(config,j);
        pool_sum += y;
        linear += y * eta[j];
      }
      
      double logprob;
      
      if(pool_sum > 0)
        logprob = std::log(1.0 - Se);
      else
        logprob = std::log(Sp);
      
      p_sub[config] = std::exp(logprob + linear);
    }
  
  //-------------------------------------------------------
    // denominator
  //-------------------------------------------------------
    
    double p_denom = std::accumulate(p_sub.begin(),
                                     p_sub.end(), 0.0);
  
  //-------------------------------------------------------
    // Calculate p_cal
  //-------------------------------------------------------
    
    NumericVector p_cal(pool_size);
  
  for(int ind = 0; ind < pool_size; ind++) {
    
    double numer = 0.0;
    
    for(int config = 0; config < n_config; config++) {
      if(y_mat(config,ind) == 1)
        numer += p_sub[config];
    }
    
    p_cal[ind] = numer / p_denom;
  }
  
  return p_cal;
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector VB_ytilde1_cpp(const NumericMatrix& x_sub,
                             const NumericVector& mu_beta,
                             const NumericMatrix& y_mat,
                             const double Se,
                             const double Sp) {
  
  int pool_size = x_sub.nrow();
  int n_config = y_mat.nrow();
  
  //-------------------------------------------------------
  // Compute x_sub %*% mu_beta
  //-------------------------------------------------------
  
  NumericVector eta(pool_size);
  
  for(int i = 0; i < pool_size; i++) {
    double temp = 0.0;
    for(int j = 0; j < x_sub.ncol(); j++)
      temp += x_sub(i,j) * mu_beta[j];
    eta[i] = temp;
  }
  
  //-------------------------------------------------------
  // Calculate p_sub
  //-------------------------------------------------------
  
  NumericVector p_sub(n_config);
  
  for(int config = 0; config < n_config; config++) {
    
    int pool_sum = 0;
    double linear = 0.0;
    
    for(int j = 0; j < pool_size; j++) {
      double y = y_mat(config,j);
      pool_sum += y;
      linear += y * eta[j];
    }
    
    double logprob;
    
    if(pool_sum > 0)
      logprob = std::log(Se);
    else
      logprob = std::log(1.0 - Sp);
    
    p_sub[config] = std::exp(logprob + linear);
  }
  
  //-------------------------------------------------------
  // Denominator
  //-------------------------------------------------------
  
  double p_denom = std::accumulate(p_sub.begin(),
                                   p_sub.end(), 0.0);
  
  //-------------------------------------------------------
  // Compute posterior probabilities
  //-------------------------------------------------------
  
  NumericVector p_cal(pool_size);
  
  for(int ind = 0; ind < pool_size; ind++) {
    
    double numer = 0.0;
    
    for(int config = 0; config < n_config; config++) {
      if(y_mat(config, ind) == 1)
        numer += p_sub[config];
    }
    
    p_cal[ind] = numer / p_denom;
  }
  
  return p_cal;
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cal_cov_sums_mpt_cpp(
    const NumericVector& mean_beta,
    const NumericVector& mean_y,
    const NumericMatrix& input_x,
    const IntegerVector& input_z,
    const double Se,
    const double Sp,
    const IntegerMatrix& pool_index,
    const NumericMatrix& y_mat){
  
  int J = pool_index.nrow();
  int pool_size = pool_index.ncol();
  int p = input_x.ncol();
  int n_config = y_mat.nrow();
  
  double sum11=0.0;
  double sum22=0.0;
  double sum33=0.0;
  double sum12=0.0;
  double sum13=0.0;
  double sum23=0.0;
  
  NumericVector eta(pool_size);
  NumericVector p_sub(n_config);
  
  for(int j=0;j<J;j++){
    
    //--------------------------------------------------
    // x_sub %*% beta
    //--------------------------------------------------
    
    for(int i=0;i<pool_size;i++){
      
      int id = pool_index(j,i)-1;      // R -> C indexing
      
      eta[i]=0.0;
      
      for(int k=0;k<p;k++)
        eta[i]+=input_x(id,k)*mean_beta[k];
    }
    
    //--------------------------------------------------
    // Compute p_sub
    //--------------------------------------------------
    
    double denom=0.0;
    
    for(int cfg=0;cfg<n_config;cfg++){
      
      int pool_sum=0;
      double linear=0.0;
      
      for(int i=0;i<pool_size;i++){
        
        double y=y_mat(cfg,i);
        
        pool_sum += y;
        linear += y*eta[i];
      }
      
      double lp;
      
      if(input_z[j]==0)
        lp=(pool_sum>0)?std::log(1-Se):std::log(Sp);
      else
        lp=(pool_sum>0)?std::log(Se):std::log(1-Sp);
      
      p_sub[cfg]=std::exp(lp+linear);
      
      denom+=p_sub[cfg];
    }
    
    //--------------------------------------------------
    // pairwise expectations
    //--------------------------------------------------
    
    for(int ii=0;ii<pool_size-1;ii++){
      
      int id1=pool_index(j,ii)-1;
      
      for(int jj=ii+1;jj<pool_size;jj++){
        
        int id2=pool_index(j,jj)-1;
        
        double Eyy=0.0;
        
        for(int cfg=0;cfg<n_config;cfg++)
          if(y_mat(cfg,ii)==1 && y_mat(cfg,jj)==1)
            Eyy+=p_sub[cfg];
          
          Eyy/=denom;
          
          double cov=
            Eyy-
            mean_y[id1]*mean_y[id2];
          
          sum11 +=
            (2.0*input_x(id1,0)*input_x(id2,0))*cov;
          
          sum22 +=
            (2.0*input_x(id1,1)*input_x(id2,1))*cov;
          
          sum33 +=
            (2.0*input_x(id1,2)*input_x(id2,2))*cov;
          
          sum12 +=
            (input_x(id1,1)+input_x(id2,1))*cov;
          
          sum13 +=
            (input_x(id1,2)+input_x(id2,2))*cov;
          
          sum23 +=
            ( input_x(id1,1)*input_x(id2,2)
                +input_x(id2,1)*input_x(id1,2) )*cov;
          
      }
    }
  }
  
  return NumericVector::create(
    sum11,
    sum22,
    sum33,
    sum12,
    sum13,
    sum23
  );
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cov_Y_mpt_cpp(
    const NumericVector& mean_beta,
    const NumericVector& mean_y,
    const NumericMatrix& input_x,
    const IntegerVector& input_z,
    const double Se,
    const double Sp,
    const IntegerMatrix& pool_index,
    const NumericMatrix& y_mat) {
  
  int J = pool_index.nrow();
  int pool_size = pool_index.ncol();
  int p = input_x.ncol();
  int n_config = y_mat.nrow();
  int n_pairs = pool_size * (pool_size - 1) / 2;
  
  NumericMatrix Eyy_cov(J, n_pairs);
  
  NumericVector eta(pool_size);
  NumericVector p_sub(n_config);
  
  for (int j = 0; j < J; j++) {
    
    //--------------------------------------------------
    // Compute x_sub %*% mean_beta
    //--------------------------------------------------
    
    for (int i = 0; i < pool_size; i++) {
      
      int id = pool_index(j, i) - 1;   // R -> C indexing
      
      eta[i] = 0.0;
      
      for (int k = 0; k < p; k++)
        eta[i] += input_x(id, k) * mean_beta[k];
    }
    
    //--------------------------------------------------
    // Compute p_sub
    //--------------------------------------------------
    
    double denom = 0.0;
    
    for (int cfg = 0; cfg < n_config; cfg++) {
      
      int pool_sum = 0;
      double linear = 0.0;
      
      for (int i = 0; i < pool_size; i++) {
        
        double y = y_mat(cfg, i);
        
        pool_sum += y;
        linear += y * eta[i];
      }
      
      double lp;
      
      if (input_z[j] == 0)
        lp = (pool_sum > 0) ? std::log(1.0 - Se) : std::log(Sp);
      else
        lp = (pool_sum > 0) ? std::log(Se) : std::log(1.0 - Sp);
      
      p_sub[cfg] = std::exp(lp + linear);
      
      denom += p_sub[cfg];
    }
    
    //--------------------------------------------------
    // Compute covariances
    //--------------------------------------------------
    
    int ct = 0;
    
    for (int ii = 0; ii < pool_size - 1; ii++) {
      
      int id1 = pool_index(j, ii) - 1;
      
      for (int jj = ii + 1; jj < pool_size; jj++) {
        
        int id2 = pool_index(j, jj) - 1;
        
        double Eyy = 0.0;
        
        for (int cfg = 0; cfg < n_config; cfg++) {
          
          if (y_mat(cfg, ii) == 1 && y_mat(cfg, jj) == 1)
            Eyy += p_sub[cfg];
        }
        
        Eyy /= denom;
        
        Eyy_cov(j, ct) = Eyy - mean_y[id1] * mean_y[id2];
        
        ct++;
      }
    }
  }
  
  return Eyy_cov;
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