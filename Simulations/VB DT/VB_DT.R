library("parallel")
library("groupTesting")
#set.seed(1)

#setwd('C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Tick Surveillance R01/Statistical Methods/Variational Bayes/')
Rcpp::sourceCpp('VB_DT.cpp')

# Compute the log-determinant of a matrix
ldet <- function(X) {
  if(!is.matrix(X)) return(log(X))
  determinant(X,logarithm = TRUE)$modulus
}

########## VB - function: Beta ###############
VB_beta <- function(input_x, input_y, input_omega, input_prior){
  PSigma_inv <- input_prior$Sigma_inv 
  Pmu_Sigma <- input_prior$mu_Sigma
  
  var_beta_inv <- crossprod(input_x*input_omega, input_x) + PSigma_inv
  var_beta <- solve(var_beta_inv)
  m_ome <- crossprod(input_x, input_y-1/2)
  mu_beta <- var_beta %*% (m_ome + Pmu_Sigma)
  return(cbind(mu_beta,var_beta))
}

############ VB_Se Function ######################
VB_Se <- function(data_pool, tilde_pool, Prior_alpha, Prior_beta){
  Post_alpha_Se <- sum(data_pool * tilde_pool) + Prior_alpha
  Post_beta_Se <- sum((1-data_pool) * tilde_pool) + Prior_beta
  Se_output <- (Post_alpha_Se)/(Post_alpha_Se + Post_beta_Se) 
  return(c(Se_output, Post_alpha_Se, Post_beta_Se))
}

############ VB_Sp Function ######################
VB_Sp <- function(data_pool, tilde_pool, Prior_alpha, Prior_beta){
  Post_alpha_Sp <- sum((1-data_pool) * (1-tilde_pool)) + Prior_alpha
  Post_beta_Sp <- sum(data_pool * (1-tilde_pool)) + Prior_beta
  Sp_output <- (Post_alpha_Sp)/(Post_alpha_Sp + Post_beta_Sp) 
  return(c(Sp_output, Post_alpha_Sp, Post_beta_Sp))
}

############# Calculate Beta Density Function #####
cal_prob_beta <- function(alpha_value, beta_value){
  gamma_cons <- lbeta(alpha_value, beta_value)
  log_beta <- digamma(alpha_value) - digamma(alpha_value + beta_value)
  log_1minbeta <- digamma(beta_value) - digamma(alpha_value + beta_value)
  kernel <- (alpha_value-1)*log_beta + (beta_value-1)*log_1minbeta
  return(-gamma_cons + kernel)
}

############ Calculate Y_tilde Joint Density Function #####
cal_density_Y0 <- function(Ey, x_abc, input_beta, 
                           Se1_alpha, Se1_beta, Sp1_alpha, Sp1_beta){
  pool_abc <- 1-prod(1-Ey)
  log_1minSe1 <- digamma(Se1_beta) - digamma(Se1_alpha+Se1_beta)
  log_Sp1 <- digamma(Sp1_alpha) - digamma(Sp1_alpha+Sp1_beta)
  
  p1 <- (pool_abc)*log_1minSe1 + (1-pool_abc)*log_Sp1 + 
    sum(Ey*x_abc%*%input_beta)
  return(p1)
}

cal_density_Y1 <- function(Ey, y_abc, x_abc,input_beta, 
                           Se1_alpha, Se1_beta, Se2_alpha, Se2_beta,
                           Sp1_alpha, Sp1_beta, Sp2_alpha, Sp2_beta){
  pool_abc <- 1-prod(1-Ey)
  log_Se1 <- digamma(Se1_alpha) - digamma(Se1_alpha + Se1_beta)
  log_1minSp1 <- digamma(Sp1_beta) - digamma(Sp1_alpha + Sp1_beta)
  log_Se2 <- digamma(Se2_alpha) - digamma(Se2_alpha + Se2_beta)
  log_1minSe2 <- digamma(Se2_beta) - digamma(Se2_alpha + Se2_beta)
  log_Sp2 <- digamma(Sp2_alpha) - digamma(Sp2_alpha + Sp2_beta)
  log_1minSp2 <- digamma(Sp2_beta) - digamma(Sp2_alpha + Sp2_beta)
  
  p1 <- pool_abc*log_Se1 + (1-pool_abc)*log_1minSp1 + 
    sum(Ey*(y_abc*log_Se2 + (1-y_abc)*log_1minSe2))+
    sum((1-Ey)*((1-y_abc)*log_Sp2 + y_abc*log_1minSp2))+
    sum(Ey*(x_abc%*%input_beta))
  return(p1)
}

product_abc <- function(a, b, c) {
  return(a * b * c)
}

#function: 1/(x*(1-x))
x_x1mins <- function(x){
  return(1/(x * (1-x)))
}

cal_Hessian_unknw <- function(mean_beta, mean_y, 
                              input_x, input_z, input_retest_y,
                              input_Se, input_Sp,
                              retest_ind,mean_z){
  p_est <- plogis(input_x %*% mean_beta)
  pq_est <- p_est*(1-p_est)
  
  J <- length(input_z)
  Ez_1mins <- 1-mean_z
  Ez <- 1-Ez_1mins # E(Z_tilde=1)
  
  Ey_sub <- mean_y[retest_ind] # E(y_tilde) from 1 to J* (in positive pools)
  Ey_sub_1mins <- 1-Ey_sub
  
  retest_y_sub <- input_retest_y[retest_ind]
  Ez_sub <- Ez[which(input_z>0)] # E(z_tilde) for positive pools only
  Ez_sub_1mins <- 1-Ez_sub
  
  p1 <- sum(-pq_est)
  p2 <- sum(-input_x[,2]^2*pq_est)
  p3 <- sum(-input_x[,3]^2*pq_est)
  
  p4 <- x_x1mins(input_Se[1])^2 * sum(-Ez*(input_Se[1]^2+input_z-2*input_z*input_Se[1])) #d2_Se_p
  p5 <- x_x1mins(input_Se[2])^2 * sum(-Ey_sub*(input_Se[2]^2+retest_y_sub-2*retest_y_sub*input_Se[2])) #d2_Se_i
  p6 <- x_x1mins(input_Sp[1])^2 * sum(-Ez_1mins*(input_Sp[1]^2+(1-input_z)-2*(1-input_z)*input_Sp[1])) #d2_Sp_p
  p7 <- x_x1mins(input_Sp[2])^2 * sum(-Ey_sub_1mins*(input_Sp[2]^2+(1-retest_y_sub)-2*(1-retest_y_sub)*input_Sp[2])) #d2_Sp_i
  
  p12 <- sum(-input_x[,2]*pq_est)
  p13 <- sum(-input_x[,3]*pq_est)
  p23 <- sum(-(input_x[,2]*input_x[,3])*pq_est)
  
  Hessian_manu <- matrix(NA, nrow=dim(input_x)+length(input_Se)+length(input_Sp), 
                         ncol=dim(input_x)+length(input_Se)+length(input_Sp))
  Hessian_manu <- diag(c(p1,p2,p3,p4,p5,p6,p7))
  Hessian_manu[1,2] <- Hessian_manu[2,1] <- p12
  Hessian_manu[1,3] <- Hessian_manu[3,1] <- p13
  Hessian_manu[2,3] <- Hessian_manu[3,2] <- p23
  
  return(Hessian_manu)
}

##### Printing current parameters ##### 
N <- 300 #Sample size

param <- c(-4,2,1) #Beta_true: prev=10%


######### Part I. Generate Data ###############
######   (Use 'groupTesting' package)   ######

## Dorfman Pool Testing (DT) - unknown Se&Sp
S <- 2 #2-stage hierarchical testing
psz <- c(10,1) #Pool sizes used in stages 1-2
Se_true <- c(0.95,0.98) #Sensitivity in stages 1-2
Sp_true <- c(0.98,0.99) #Specificity in stages 1-2
assayID <- c(1,1) #Assays used in stages 1-2
tau_sq = 100

y.mat <- as.matrix(expand.grid(rep(list(0:1),psz[1])))
row.index.mat <- matrix(NA,dim(y.mat)[2],dim(y.mat)[1])
for(r in 1:dim(y.mat)[2]){
  row.index.mat[r,] = y.mat[,r] == 1
}

# Define the indices for combinations
indices <- t(combn(c(1:length(param)),2))

con_tol_vec <- c(1,0.1,0.01,10^-3,10^-4,10^-5,10^-6,10^-7,10^-8) #converge criteria
con_tol <- con_tol_vec[1]
maxiter <- 10000 #maximum iterations in VB

n_sim <- 1

simu_func <- function(seed){

  results.list <- list()
  con.ct <- 1
  
  x <- cbind(1, rnorm(N), rbinom(N,1,0.5))
  colnames(x) <- c("Intercept", "x1", "x2")
  
  pReg <- plogis(x %*% param)
  # Simulating test responses
  DT.gtData <- hier.gt.simulation(N, p=pReg, S, psz, Se=Se_true, Sp=Sp_true, assayID)$gtData
  
  ############## Part II. VB Data Preparation ############
  size <- psz[1] # Pool size
  J <- N/size #number of pools
  pool_index_mtx <- DT.gtData[c(1:J),c(6:(6+psz[1]-1))]
  data_x <- x
  data_z <- DT.gtData[c(1:J),1]
  retest_y <- rep(0,N) #initialize retesting info - all 0s
  retest_ind <- DT.gtData[c((J+1):dim(DT.gtData)[1]),6] #index of retesting individuals
  retest_y[retest_ind] <- DT.gtData[c((J+1):dim(DT.gtData)[1]),1]
  retest_y_sub <- DT.gtData[c((J+1):dim(DT.gtData)[1]),1] #retesting info only
  p <- dim(data_x)[2] # p
  
  ############ CAVI algorithm ################
  time1 <- Sys.time()

  lowerbound <- numeric(maxiter)
  PSigma <- tau_sq*diag(p) #Prior of Sigma_beta
  PSigma_inv <- solve(PSigma) 
  Pmu <- rep(0,p) #Prior of mu_beta
  Pmu_Sigma <- c(PSigma_inv %*% Pmu)
  Pri_beta <- list(Sigma_inv = PSigma_inv, mu_Sigma=Pmu_Sigma)
  Pdet <- ldet(PSigma_inv) 
  Palpha <- 1 #Prior of Se&Sp: alpha(e), alpha(p) 
  Pbeta <- 1 #Prior of Se&Sp: beta(e), beta(p)
  
  # Initialization
  omega0 <- rep(1/4, N)
  y0 <- rep(0, N)
  alpha0 <- c(rnorm(1,703,10),rnorm(1,763,10),rnorm(1,3197,10),rnorm(1,1109,10)) #alpha0 for Se_p,Se_i,Sp_p,Sp_i
  beta0 <- c(rnorm(1,39,10),rnorm(1,22,10),rnorm(1,65,10),rnorm(1,13,10)) #beta0 for Se_p,Se_i,Sp_p,Sp_i
  beta0[beta0<0] <- 0.1
  Se0_1 <- c(alpha0[1]/(alpha0[1] + beta0[1]), alpha0[1], beta0[1])
  Sp0_1 <- c(alpha0[3]/(alpha0[3] + beta0[3]), alpha0[3], beta0[3])
  Se0_2 <- c(alpha0[2]/(alpha0[2] + beta0[2]), alpha0[2], beta0[2])
  Sp0_2 <- c(alpha0[4]/(alpha0[4] + beta0[4]), alpha0[4], beta0[4])
  # t=1
  # Update Beta
  beta1 <- VB_beta(data_x, y0, omega0, Pri_beta) #output: mu_beta + var_beta
  mu_beta <- c(beta1[,1])
  var_beta <- beta1[,2:(p+1)]
  # Update Omega
  eta <- c(data_x %*% mu_beta)
  xi <- sqrt(eta^2 +  rowSums(data_x %*% var_beta * data_x))
  omega1 <- tanh(xi/2)/(2*xi)
  # Update y_tilde & Z_tilde
  vq_z <- rep(0, J)
  pool_data <- rep(NA, size)
  z0_a <- c() #Z_tilde for pool response
  z0_b <- c() #Z_tilde for the individuals in positive pools
  for (j in 1:J){
    pool_index <- pool_index_mtx[j,]
    x_sub <- data_x[pool_index,]
    if (data_z[j] == 0){
      pool_data <- VB_ytilde0_cpp(x_sub, mu_beta,y.mat, Se0_1[2], Se0_1[3], Sp0_1[2], Sp0_1[3])
      vq_z[j] <- cal_density_Y0(pool_data[1:psz[1]], 
                                x_sub, mu_beta, Se0_1[2], Se0_1[3], Sp0_1[2], Sp0_1[3])
      y0[pool_index] <- pool_data[1:psz[1]]
      z0_a <- c(z0_a, 1-pool_data[(psz[1]+1)])
    }
    if (data_z[j] == 1){
      y_sub <- retest_y[pool_index]
      pool_data <- VB_ytilde1_cpp(y_sub, x_sub, mu_beta,y.mat, Se0_1[2], Se0_1[3], Se0_2[2], Se0_2[3], Sp0_1[2], Sp0_1[3], Sp0_2[2], Sp0_2[3])
      vq_z[j]  <- cal_density_Y1(pool_data[1:psz[1]], 
                           y_sub, x_sub, mu_beta, Se0_1[2], Se0_1[3], Se0_2[2], Se0_2[3], Sp0_1[2], Sp0_1[3], Sp0_2[2], Sp0_2[3])
      y0[pool_index] <- pool_data[1:psz[1]]
      z0_a <- c(z0_a, 1-pool_data[(psz[1]+1)])
      z0_b <- c(z0_b, pool_data[1:psz[1]])
    }
  }
  #Update Se&Spe
  Se0_1 <- VB_Se(data_z, z0_a, Palpha, Pbeta) #pool sensitivity
  Se0_2 <- VB_Se(retest_y_sub, z0_b, Palpha, Pbeta) #individual pool sensitivity
  Sp0_1 <- VB_Sp(data_z, z0_a, Palpha, Pbeta) #pool specificity
  Sp0_2 <- VB_Sp(retest_y_sub, z0_b, Palpha, Pbeta) #individual pool sensitivity
  
  lowerbound[1] <- sum(vq_z) + 
    cal_prob_beta(Se0_1[2], Se0_1[3]) +
    cal_prob_beta(Se0_2[2], Se0_2[3]) +
    cal_prob_beta(Sp0_1[2], Sp0_1[3]) +
    cal_prob_beta(Sp0_2[2], Sp0_2[3]) +
    0.5*p + 0.5*ldet(var_beta) + 0.5*Pdet - 0.5*t(mu_beta - Pmu)%*%PSigma_inv%*%(mu_beta - Pmu) + 
    sum((y0-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(PSigma_inv %*% var_beta))
  
  # Iteration t=2:maxiter
  output.VB <- list()#Beta(mu+var)+Se+Sp+Ey
  for(t in 2:maxiter){
    
    # Update of beta
    beta1 <- VB_beta(data_x, y0, omega1, Pri_beta) #output: mu_beta + var_beta
    mu_beta <- c(beta1[,1])
    var_beta <- beta1[,2:(p+1)]
    
    # Update of omega
    eta <- c(data_x %*% mu_beta)
    xi <- sqrt(eta^2 +  rowSums(data_x %*% var_beta * data_x))
    omega1 <- tanh(xi/2)/(2*xi)
    
    # Update y_tilde & Calculate q(Z_tilde)
    vq_z <- rep(0, J)
    z0_a <- c() #Z_tilde for pool response
    z0_b <- c() #Z_tilde for the individuals in positive pools
    for (j in 1:J){
      pool_index <- pool_index_mtx[j,]
      x_sub <- data_x[pool_index,]
      if (data_z[j] == 0){
        pool_data <- VB_ytilde0_cpp(x_sub, mu_beta,y.mat, Se0_1[2], Se0_1[3], Sp0_1[2], Sp0_1[3])
        vq_z[j] <- cal_density_Y0(pool_data[1:psz[1]], 
                                  x_sub, mu_beta,Se0_1[2], Se0_1[3], Sp0_1[2], Sp0_1[3])
        y0[pool_index] <- pool_data[1:psz[1]]
        z0_a <- c(z0_a, 1-pool_data[(psz[1]+1)])
      }
      if (data_z[j] == 1){
        y_sub <- retest_y[pool_index]
        pool_data <- VB_ytilde1_cpp(y_sub, x_sub, mu_beta,y.mat, Se0_1[2], Se0_1[3], Se0_2[2], Se0_2[3], Sp0_1[2], Sp0_1[3], Sp0_2[2], Sp0_2[3])
        vq_z[j] <- cal_density_Y1(pool_data[1:psz[1]],
                                  y_sub, x_sub, mu_beta, Se0_1[2], Se0_1[3], Se0_2[2], Se0_2[3], Sp0_1[2], Sp0_1[3], Sp0_2[2], Sp0_2[3])
        y0[pool_index] <- pool_data[1:psz[1]]
        z0_a <- c(z0_a, 1-pool_data[(psz[1]+1)])
        z0_b <- c(z0_b, pool_data[1:psz[1]])
      }
    }
    
    #Update Se&Spe
    Se0_1 <- VB_Se(data_z, z0_a, Palpha, Pbeta) #pool sensitivity
    Se0_2 <- VB_Se(retest_y_sub, z0_b, Palpha, Pbeta) #individual pool sensitivity
    Sp0_1 <- VB_Sp(data_z, z0_a, Palpha, Pbeta) #pool specificity
    Sp0_2 <- VB_Sp(retest_y_sub, z0_b, Palpha, Pbeta) #individual pool sensitivity
    
    lowerbound[t]  <- sum(vq_z) + 
      cal_prob_beta(Se0_1[2], Se0_1[3]) +
      cal_prob_beta(Se0_2[2], Se0_2[3]) +
      cal_prob_beta(Sp0_1[2], Sp0_1[3]) +
      cal_prob_beta(Sp0_2[2], Sp0_2[3]) +
      0.5*p + 0.5*ldet(var_beta) + 0.5*Pdet - 0.5*t(mu_beta - Pmu)%*%PSigma_inv%*%(mu_beta - Pmu) + 
      sum((y0-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(PSigma_inv %*% var_beta))
    
    if(abs(lowerbound[t] - lowerbound[t-1]) < con_tol) {
      results <- list(mu_beta=mu_beta,
                      var_beta=beta1[,2:(1+p)],
                      Se=c(Se0_1[1], Se0_2[1]),
                      Sp=c(Sp0_1[1], Sp0_2[1]),
                      Se_alp_beta=c(Se0_1[2], Se0_1[3], Se0_2[2], Se0_2[3]), #Se_p:alpha, beta, Se_i: alpha, beta
                      Sp_alp_beta=c(Sp0_1[2], Sp0_1[3], Sp0_2[2], Sp0_2[3]), #Sp_p:alpha, beta, Sp_i: alpha, beta
                      est_y=y0,est_z = z0_a)
      
      Ez_1mins <- 1-results$est_z
      Ez <- 1-Ez_1mins # E(Z_tilde=1)
      
      Ey_sub <- results$est_y[retest_ind] # E(y_tilde) from 1 to J* (in positive pools)
      Ey_sub_1mins <- 1-Ey_sub
      
      Ez_sub <- Ez[which(data_z>0)] # E(z_tilde) for positive pools only
      Ez_sub_1mins <- 1-Ez_sub
      Hessian_manu <- cal_Hessian_unknw(results$mu_beta, results$est_y, 
                                        data_x, data_z, retest_y,
                                        results$Se, results$Sp,
                                        retest_ind, results$est_z)
      
      Ey_Ey1mins <- results$est_y * (1-results$est_y)
      cov_sums <- cal_cov_sums_unknw_cpp(results$mu_beta, results$est_y, 
                                     data_x, data_z, retest_y, 
                                     results$Se, results$Sp, 
                                     results$Se_alp_beta, results$Sp_alp_beta, 
                                     pool_index_mtx, results$est_z,y.mat)
      var_dbeta <- colSums(data_x^2*Ey_Ey1mins) + cov_sums[1:p]
      
      var_Se1 <- (x_x1mins(results$Se[1])^2) * sum((data_z-results$Se[1])^2 * Ez * Ez_1mins)
      var_Sp1 <- (x_x1mins(results$Sp[1])^2) * sum((1-data_z-results$Sp[1])^2 * Ez * Ez_1mins)
      var_Se2 <- (x_x1mins(results$Se[2])^2) * (sum((retest_y_sub-results$Se[2])^2*Ey_sub*Ey_sub_1mins) +
                                                  cov_sums[p+1])
      var_Sp2 <- (x_x1mins(results$Sp[2])^2) * (sum((1-retest_y_sub-results$Sp[2])^2*Ey_sub*Ey_sub_1mins) +
                                                  cov_sums[p+2])
      
      cov_par <- diag(c(var_dbeta, var_Se1, var_Se2, var_Sp1, var_Sp2))
      
      cov_d0d1 <- sum(data_x[,1]*data_x[,2]*Ey_Ey1mins) + cov_sums[p+3] #d1: beta0 vs beta1
      cov_d0d2 <- sum(data_x[,1]*data_x[,3]*Ey_Ey1mins) + cov_sums[p+4] #d1: beta0 vs beta2
      cov_d1d2 <- sum(data_x[,2]*data_x[,3]*Ey_Ey1mins) + cov_sums[p+5] #d1: beta1 vs beta2
      cov_dbeta <- c(cov_d0d1, cov_d0d2, cov_d1d2)
      
      cov_Se1_Se2 <- x_x1mins(results$Se[1]) * x_x1mins(results$Se[2]) *
        sum((1-results$Se[1]) * (retest_y_sub-results$Se[2]) * 
              Ey_sub * rep(Ez_sub_1mins,each=size))
      
      cov_Se1_Sp1 <- x_x1mins(results$Se[1]) * x_x1mins(results$Sp[1]) * 
        sum((data_z-results$Se[1])*(1-data_z-results$Sp[1]) * (-Ez*Ez_1mins))
      
      cov_Se1_Sp2 <- x_x1mins(results$Se[1]) * x_x1mins(results$Sp[2]) *
        sum((1-results$Se[1]) * (1-retest_y_sub-results$Sp[2]) *
              Ey_sub * rep(-Ez_sub_1mins, each=size))
      
      cov_Se2_Sp1 <- x_x1mins(results$Se[2]) * x_x1mins(results$Sp[1]) *
        sum((0-results$Sp[1]) * (retest_y_sub-results$Se[2]) *
              rep(-Ez_sub_1mins, each=size) * Ey_sub)
      
      cov_Sp1_Sp2 <- x_x1mins(results$Sp[1]) * x_x1mins(results$Sp[2]) *
        sum((-results$Sp[1]) * (1-retest_y_sub-results$Sp[2])* 
              rep(Ez_sub_1mins,each=size) * Ey_sub)
      
      cov_Se1_par <- x_x1mins(results$Se[1]) * cov_sums[(p+6):(p+5+p)] #d_Se_p vs Beta's
      cov_Sp1_par <- x_x1mins(results$Sp[1]) * cov_sums[(p+6+p):(p+5+2*p)] #d_Sp_p vs Beta's
      cov_Se2_par <- x_x1mins(results$Se[2]) * cov_sums[(p+6+2*p):(p+5+3*p)] #d_Se_i vs Beta's
      cov_Sp2_par <- x_x1mins(results$Sp[2]) * cov_sums[(p+6+3*p):(p+5+4*p)] #d_Sp_i vs Beta's
      
      cov_Se2_Sp2 <- x_x1mins(results$Se[2]) * x_x1mins(results$Sp[2]) *  cov_sums[p+6+4*p] #d_Se_i vs d_Sp_i
      
      # Assign cov_xy to cov_par(beta part) using loop
      for (i in 1:nrow(indices)) {
        cov_par[indices[i, 1], indices[i, 2]] <- cov_par[indices[i, 2], indices[i, 1]] <- cov_dbeta[i]
      }
      # Assign cov_SeSp to cov_par(Se&Sp part) manually
      cov_par[4,5] <- cov_par[5,4] <- cov_Se1_Se2
      cov_par[4,6] <- cov_par[6,4] <- cov_Se1_Sp1
      cov_par[4,7] <- cov_par[7,4] <- cov_Se1_Sp2
      
      cov_par[5,6] <- cov_par[6,5] <- cov_Se2_Sp1
      cov_par[5,7] <- cov_par[7,5] <- cov_Se2_Sp2
      
      cov_par[6,7] <- cov_par[7,6] <- cov_Sp1_Sp2
      
      # Assign cov_SeSp_par manually
      cov_par[4,c(1:p)] <- cov_par[c(1:p),4] <- cov_Se1_par
      cov_par[5,c(1:p)] <- cov_par[c(1:p),5] <- cov_Se2_par
      cov_par[6,c(1:p)] <- cov_par[c(1:p),6] <- cov_Sp1_par
      cov_par[7,c(1:p)] <- cov_par[c(1:p),7] <- cov_Sp2_par
      
      cov_VB_par <- solve(-Hessian_manu - cov_par)
      
      
      time2 <- Sys.time()
      time.diff <- as.numeric(difftime(time2, time1, units = "secs"))
      
      VB.mean <- c(results$mu_beta, results$Se, results$Sp)
      #if (sum(is.na(sqrt(diag(cov_VB_par)))) > 0) {
        #print(paste0('Seed ',seed, ": NaNs produced"))}
      VB.var <- c(diag(results$var_beta),
                  (results$Se_alp_beta[1]*results$Se_alp_beta[2])/((results$Se_alp_beta[1]+results$Se_alp_beta[2])^2*(results$Se_alp_beta[1]+results$Se_alp_beta[2]+1)),
                  (results$Se_alp_beta[3]*results$Se_alp_beta[4])/((results$Se_alp_beta[3]+results$Se_alp_beta[4])^2*(results$Se_alp_beta[3]+results$Se_alp_beta[4]+1)),
                  (results$Sp_alp_beta[1]*results$Sp_alp_beta[2])/((results$Sp_alp_beta[1]+results$Sp_alp_beta[2])^2*(results$Sp_alp_beta[1]+results$Sp_alp_beta[2]+1)),
                  (results$Sp_alp_beta[3]*results$Sp_alp_beta[4])/((results$Sp_alp_beta[3]+results$Sp_alp_beta[4])^2*(results$Sp_alp_beta[3]+results$Sp_alp_beta[4]+1)))
      VB.std <- sqrt(VB.var)
      VB.lower <- VB.mean - 1.96*VB.std
      VB.upper <- VB.mean + 1.96*VB.std
      VB.coverage <- VB.lower <= c(param,Se_true,Sp_true) & VB.upper >= c(param,Se_true,Sp_true)
      
      VB.std_cor <- sqrt(abs(diag(cov_VB_par))) 
      VB.lower_cor <- VB.mean - 1.96*VB.std_cor
      VB.upper_cor <- VB.mean + 1.96*VB.std_cor
      VB.coverage_cor <- VB.lower_cor <= c(param,Se_true,Sp_true) & 
      VB.upper_cor >= c(param,Se_true,Sp_true)
      
      VB.alt.low <- c(0,0,0,
                      qbeta(0.025,results$Se_alp_beta[1],results$Se_alp_beta[2]),
                      qbeta(0.025,results$Se_alp_beta[3],results$Se_alp_beta[4]),
                      qbeta(0.025,results$Sp_alp_beta[1],results$Sp_alp_beta[2]),
                      qbeta(0.025,results$Sp_alp_beta[3],results$Sp_alp_beta[4])
                      )
      VB.alt.upp <- c(0,0,0,
                      qbeta(0.975,results$Se_alp_beta[1],results$Se_alp_beta[2]),
                      qbeta(0.975,results$Se_alp_beta[3],results$Se_alp_beta[4]),
                      qbeta(0.975,results$Sp_alp_beta[1],results$Sp_alp_beta[2]),
                      qbeta(0.975,results$Sp_alp_beta[3],results$Sp_alp_beta[4])
      )
      VB.alt.cov <- VB.alt.low <= c(param,Se_true,Sp_true) & VB.alt.upp >= c(param,Se_true,Sp_true)
      
      VB.output <- data.frame(true_pars = c(param,Se_true,Sp_true),
                              est_mean = VB.mean,
                              est_std = VB.std,
                              est_lci = VB.lower,
                              est_uci = VB.upper,
                              coverage = VB.coverage,
                              est_std_cor = VB.std_cor,
                              est_lci_cor = VB.lower_cor,
                              est_uci_cor = VB.upper_cor,
                              coverage_cor = VB.coverage_cor,
                              est_lci_alt = VB.alt.low, 
                              est_uci_alt = VB.alt.upp,
                              est_alt_cov = VB.alt.cov,
                              time = time.diff)
      results.list[[con.ct]] <- VB.output
      print(con_tol)
      if(con.ct < length(con_tol_vec)){
        con.ct = con.ct + 1
        con_tol <- con_tol_vec[con.ct]
      }else{
        break
      }
    }
    
  }
  return(results.list)
}

out.mat <- NULL
time.vec <- NULL
for(sims in 1:n_sim){
  output.cmb = simu_func(1)
  for(j in 1:length(output.cmb)){
    output = output.cmb[[j]]
    out.mat <- rbind(out.mat,output[[1]], output[[2]], output[[3]], output[[4]], output[[5]], output[[6]], output[[7]], output[[8]], output[[9]], output[[10]], output[[11]], output[[12]], output[[13]], output[[14]])
    time.vec <- c(time.vec,output[[11]][1])
  }
}



write.table(out.mat,'betares.txt',row.names = FALSE, col.names = FALSE) 
write.table(time.vec,'timeres.txt',row.names = FALSE, col.names = FALSE)
