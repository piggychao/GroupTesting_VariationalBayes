library("parallel")
library("groupTesting")

rm(list=ls())
source("VBFuncUnk.txt")

nm_cores <- 60

##### Capture arguments passed from command line/bash script ##### 
args <- commandArgs(trailingOnly = TRUE)

nm_set <- as.numeric(args)+18 #randomly choose a set number

if(length(args)==1){
  ##### Create settings matrix ##### 
  settings <- expand.grid(n = c(3000, 6000, 12000), 
                          prev = c(10, 34)) 
  args <- as.numeric(settings[as.numeric(args),])
}

##### Printing current parameters ##### 
N <- as.numeric(args[1]) #Sample size

if (as.numeric(args[2])==10){
  param <- c(-4,2,1) #Beta_true: prev=10%
}
if (as.numeric(args[2])==34){
  param <- c(-2,3,1) #Beta_true: prev=34%
}

######### Part I. Generate Data ###############
######   (Use 'groupTesting' package)   ######

## Dorfman Pool Testing (DT) - unknown Se&Sp
S <- 2 #2-stage hierarchical testing
psz <- c(3,1) #Pool sizes used in stages 1-2
Se_true <- c(0.95,0.98) #Sensitivity in stages 1-2
Sp_true <- c(0.98,0.99) #Specificity in stages 1-2
assayID <- c(1,1) #Assays used in stages 1-2

# Define the indices for combinations
indices <- t(combn(c(1:length(param)),2))

con_tol <- 1e-8 #converge criteria
maxiter <- 10000 #maximum iterations in VB

n_sim <- 500

simu_func <- function(seed){
  set.seed(seed)
  
  x <- cbind(1, rnorm(N), rbinom(N,1,0.5))
  colnames(x) <- c("Intercept", "x1", "x2")
  
  pReg <- plogis(x %*% param)
  # Simulating test responses
  DT.gtData <- hier.gt.simulation(N, p=pReg, S, psz, Se=Se_true, Sp=Sp_true, assayID)$gtData
  
  ############## Part II. VB Data Preparation ############
  size <- psz[1] # Pool size
  J <- N/size #number of pools
  pool_index_mtx <- DT.gtData[c(1:J),c(6:8)]
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
  PSigma <- diag(p) #Prior of Sigma_beta
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
      pool_data <- VB_ytilde0(x_sub, mu_beta, alpha0[1], beta0[1], alpha0[3], beta0[3])
      vq_z[j] <- cal_density_Y0(pool_data[1], pool_data[2], pool_data[3], 
                                x_sub, mu_beta, alpha0[1], beta0[1], alpha0[3], beta0[3])
      y0[pool_index] <- c(pool_data[1], pool_data[2], pool_data[3])
      z0_a <- c(z0_a, 1-pool_data[4])
    }
    if (data_z[j] == 1){
      y_sub <- retest_y[pool_index]
      pool_data <- VB_ytilde1(y_sub, x_sub, mu_beta, alpha0[1], beta0[1], alpha0[2], beta0[2], alpha0[3], beta0[3], alpha0[4], beta0[4])
      vq_z[j]  <- cal_density_Y1(pool_data[1], pool_data[2], pool_data[3], 
                                 y_sub, x_sub, mu_beta, alpha0[1], beta0[1], alpha0[2], beta0[2], alpha0[3], beta0[3], alpha0[4], beta0[4])
      y0[pool_index] <- c(pool_data[1], pool_data[2], pool_data[3])
      z0_a <- c(z0_a, 1-pool_data[4])
      z0_b <- c(z0_b, c(pool_data[1], pool_data[2], pool_data[3]))
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
  results <- list() #Beta(mu+var)+Se+Sp+Ey
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
        pool_data <- VB_ytilde0(x_sub, mu_beta, Se0_1[2], Se0_1[3], Sp0_1[2], Sp0_1[3])
        vq_z[j] <- cal_density_Y0(pool_data[1], pool_data[2], pool_data[3], 
                                  x_sub, mu_beta, alpha0[1], beta0[1], alpha0[3], beta0[3])
        y0[pool_index] <- c(pool_data[1], pool_data[2], pool_data[3])
        z0_a <- c(z0_a, 1-pool_data[4])
      }
      if (data_z[j] == 1){
        y_sub <- retest_y[pool_index]
        pool_data <- VB_ytilde1(y_sub, x_sub, mu_beta, Se0_1[2], Se0_1[3], Se0_2[2], Se0_2[3], Sp0_1[2], Sp0_1[3], Sp0_2[2], Sp0_2[3])
        vq_z[j] <- cal_density_Y1(pool_data[1], pool_data[2], pool_data[3], 
                                  y_sub, x_sub, mu_beta, alpha0[1], beta0[1], alpha0[2], beta0[2], alpha0[3], beta0[3], alpha0[4], beta0[4])
        y0[pool_index] <- c(pool_data[1], pool_data[2], pool_data[3])
        z0_a <- c(z0_a, 1-pool_data[4])
        z0_b <- c(z0_b, c(pool_data[1], pool_data[2], pool_data[3]))
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
                      est_y=y0)
      break
    }
    
   # print(paste0(t, ": ",lowerbound[t], " , ", abs(lowerbound[t] - lowerbound[t-1]), " , ", Se0_1[1], " , ",Se0_2[1], " , ",Sp0_1[1], " , ", Sp0_2[1]))
  }
 
  Ez.mtx <- matrix(1-results$est_y, nrow = J, byrow=T) 
  Ez_1mins <- mapply(product_abc, Ez.mtx[,1], Ez.mtx[,2], Ez.mtx[,3]) # E(Z_tilde=0)=1-E(Z_tilde=1)
  Ez <- 1-Ez_1mins # E(Z_tilde=1)
  
  Ey_sub <- results$est_y[retest_ind] # E(y_tilde) from 1 to J* (in positive pools)
  Ey_sub_1mins <- 1-Ey_sub
  
  Ez_sub <- Ez[which(data_z>0)] # E(z_tilde) for positive pools only
  Ez_sub_1mins <- 1-Ez_sub
  Hessian_manu <- cal_Hessian_unknw(results$mu_beta, results$est_y, 
                                    data_x, data_z, retest_y,
                                    results$Se, results$Sp,
                                    retest_ind)
  
  Ey_Ey1mins <- results$est_y * (1-results$est_y)
  cov_sums <- cal_cov_sums_unknw(results$mu_beta, results$est_y, 
                                 data_x, data_z, retest_y, 
                                 results$Se, results$Sp, 
                                 results$Se_alp_beta, results$Sp_alp_beta, 
                                 pool_index_mtx)
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

  cov_Se2_Sp2 <- x_x1mins(results$Se[2]) * x_x1mins(results$Sp[2]) * 
      (sum((retest_y_sub-results$Se[2])*(1-retest_y_sub-results$Sp[2])*
           (-Ey_sub* Ey_sub_1mins)) + cov_sums[p+6+4*p]) #d_Se_i vs d_Sp_i
  
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
  if (sum(is.na(sqrt(diag(cov_VB_par)))) > 0) {
    print(paste0('Seed ',seed, ": NaNs produced"))}
  VB.std <- c(sqrt(diag(results$var_beta)),
              (results$Se_alp_beta[1]*results$Se_alp_beta[2])/((results$Se_alp_beta[1]+results$Se_alp_beta[2])^2*(results$Se_alp_beta[1]+results$Se_alp_beta[2]+1)),
              (results$Se_alp_beta[3]*results$Se_alp_beta[4])/((results$Se_alp_beta[3]+results$Se_alp_beta[4])^2*(results$Se_alp_beta[3]+results$Se_alp_beta[4]+1)),
              (results$Sp_alp_beta[1]*results$Sp_alp_beta[2])/((results$Sp_alp_beta[1]+results$Sp_alp_beta[2])^2*(results$Sp_alp_beta[1]+results$Sp_alp_beta[2]+1)),
              (results$Sp_alp_beta[1]*results$Sp_alp_beta[2])/((results$Sp_alp_beta[3]+results$Sp_alp_beta[4])^2*(results$Sp_alp_beta[3]+results$Sp_alp_beta[4]+1)))
  VB.lower <- VB.mean - 1.96*VB.std
  VB.upper <- VB.mean + 1.96*VB.std
  VB.coverage <- VB.lower <= c(param,Se_true,Sp_true) & VB.upper >= c(param,Se_true,Sp_true)
  
  VB.std_cor <- sqrt(abs(diag(cov_VB_par))) 
  VB.lower_cor <- VB.mean - 1.96*VB.std_cor
  VB.upper_cor <- VB.mean + 1.96*VB.std_cor
  VB.coverage_cor <- VB.lower_cor <= c(param,Se_true,Sp_true) & 
                     VB.upper_cor >= c(param,Se_true,Sp_true)
  
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
                          seed = seed,
                          time = time.diff)
  
  return(VB.output)
}

output <- mclapply(((nm_set-1)*n_sim+1):(nm_set*n_sim), 
                   simu_func, mc.cores = nm_cores)
time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

saveRDS(output, file = paste0('VB_DT_unknown_', time,'.rds'))
  



