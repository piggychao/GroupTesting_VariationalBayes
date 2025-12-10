library("parallel")
library("groupTesting")

rm(list=ls())
source("VBFuncMPT.txt")
nm_cores <- 60

##### Capture arguments passed from command line/bash script ##### 
args <- commandArgs(trailingOnly = TRUE)

nm_set <- as.numeric(args)+6 #randomly choose a set number

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

## Master Pool Testing (MPT) - known Se&Sp
S <- 1 #1-stage, master pool testing
psz <- 3 #pool size
Se_true <- 0.95 #Pool sensitivity
Sp_true <- 0.98 #Pool specificity
assayID <- 1 #Assays used in master pools

# Define the indices for combinations
indices <- t(combn(c(1:length(param)),2))

con_tol <- 1e-8 #converge criteria
maxiter <- 10000 #maximum iterations in VB

n_sim <- 500

simu_func <- function(seed){
 set.seed(seed)
 ## Model M2  
 ## Group testing data generation
 x <- cbind(1, rnorm(N), rbinom(N,1,0.5))
 colnames(x) <- c("Intercept", "x1", "x2")

 pReg <- plogis(x %*% param)
 # Simulating test responses
 MPT.gtData <- hier.gt.simulation(N, p=pReg, S, psz, Se=Se_true, Sp=Sp_true, assayID)$gtData
 
 J <- N/psz[1]
 pool_index_mtx <- MPT.gtData[c(1:J),c(6:8)]
 data_x <- x
 data_z <- MPT.gtData[c(1:J),1]
 p <- dim(data_x)[2] # p
 size <- psz[1]      # Pool size
 lowerbound <- numeric(maxiter) #initial value of the ELBO in VB
 
 time1 <- Sys.time()
 ######### Part I. VB CAVI algorithm ################
 PSigma <- diag(p) #Prior of Sigma_beta
 PSigma_inv <- solve(PSigma) 
 Pmu <- c(rep(0,p)) #Prior of mu_beta
 Pmu_Sigma <- c(PSigma_inv %*% Pmu)
 Pbeta <- list(Sigma_inv = PSigma_inv, mu_Sigma=Pmu_Sigma)
 Pdet <- ldet(PSigma_inv) 
 
 # Initialization
 omega0 <- rep(1/4, N)
 y0 <- rep(0, N)
 
 # t=1
 # Update Beta
 beta1 <- VB_beta(data_x, y0, omega0, Pbeta) #output: mu_beta + var_beta
 mu_beta <- c(beta1[,1])
 var_beta <- beta1[,2:(p+1)]
 # Update Omega
 eta <- c(data_x %*% mu_beta)
 xi <- sqrt(eta^2 +  rowSums(data_x %*% var_beta * data_x))
 omega1 <- tanh(xi/2)/(2*xi)
 # Update y_tilde
 vq_z <- rep(0, J)
 pool_data <- rep(NA, size)
 for (j in 1:J){
   pool_index <- pool_index_mtx[j,]
   x_sub <- data_x[pool_index,]
   if (data_z[j] == 0){
     pool_data <- VB_ytilde0(x_sub, mu_beta, Se_true, Sp_true)
     vq_z[j] <- cal_density_Y0(pool_data[1], pool_data[2], pool_data[3], 
                               x_sub, mu_beta,  Se_true, Sp_true)
     y0[pool_index] <- pool_data
   }
   if (data_z[j] == 1){
     pool_data <- VB_ytilde1(x_sub, mu_beta, Se_true, Sp_true)
     vq_z[j] <- cal_density_Y1(pool_data[1], pool_data[2], pool_data[3], 
                               x_sub, mu_beta, Se_true, Sp_true)
     y0[pool_index] <- pool_data
   }
 }
 
 lowerbound[1] <- sum(vq_z) + 
   0.5*p + 0.5*ldet(var_beta) + 0.5*Pdet - 0.5*t(mu_beta - Pmu)%*%PSigma_inv%*%(mu_beta - Pmu) + 
   sum((y0-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(PSigma_inv %*% var_beta))
 
 # Iteration t=2:maxiter
 output.VB <- list()
 for(t in 2:maxiter){
   
   # Update of beta
   beta1 <- VB_beta(data_x, y0, omega1, Pbeta) #output: mu_beta + var_beta
   mu_beta <- c(beta1[,1])
   var_beta <- beta1[,2:(p+1)]
   
   # Update of omega
   eta <- c(data_x %*% mu_beta)
   xi <- sqrt(eta^2 +  rowSums(data_x %*% var_beta * data_x))
   omega1 <- tanh(xi/2)/(2*xi)
   
   # Update y_tilde & Calculate q(Z_tilde)
   vq_z <- rep(0, J)
   for (j in 1:J){
     pool_index <- pool_index_mtx[j,]
     x_sub <- data_x[pool_index,]
     if (data_z[j] == 0){
       pool_data <- VB_ytilde0(x_sub, mu_beta, Se_true, Sp_true)
       vq_z[j] <- cal_density_Y0(pool_data[1], pool_data[2], pool_data[3], 
                                 x_sub, mu_beta,  Se_true, Sp_true)
       y0[pool_index] <- pool_data
     }
     if (data_z[j] == 1){
       pool_data <- VB_ytilde1(x_sub, mu_beta, Se_true, Sp_true)
       vq_z[j] <- cal_density_Y1(pool_data[1], pool_data[2], pool_data[3], 
                                 x_sub, mu_beta, Se_true, Sp_true)
       y0[pool_index] <- pool_data
     }
   }
   
   lowerbound[t]  <- sum(vq_z) + 
     0.5*p + 0.5*ldet(var_beta) + 0.5*Pdet - 0.5*t(mu_beta - Pmu)%*%PSigma_inv%*%(mu_beta - Pmu) + 
     sum((y0-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(PSigma_inv %*% var_beta))
   
   if(abs(lowerbound[t] - lowerbound[t-1]) < con_tol) {
     output.VB <- list(mu_beta=mu_beta, 
                       var_beta=var_beta, 
                       est_y=y0)
     break
   }
  
 }
 
 ## H Manually Calculation
 Hessian_manu <- cal_Hessian_it(output.VB$mu_beta, output.VB$est_y, data_x)
 ## Cov of first derivatives_beta calculation
 Cov_sums <- cal_cov_sums_mpt(output.VB$mu_beta, output.VB$est_y, 
                              data_x, data_z, Se_true, Sp_true, pool_index_mtx)
 
 Ey_1mins <- output.VB$est_y * (1-output.VB$est_y)
 Var_dbeta <- colSums(data_x^2*Ey_1mins) + Cov_sums[1:p] #Variance of first derivatives: beta0, beta1, beta2
 # Covariance of first derivatives: beta0-1, beta0-2, beta1-2
 multipliers <- apply(indices, 1, 
                        function(vec) data_x[,vec[1]] * data_x[, vec[2]])
 Cov_dbeta <- colSums(multipliers*Ey_1mins) + Cov_sums[c(p+1):(2*p)]
 
 cov_beta <- diag(Var_dbeta)
 # Assign values from cov_xy to cov_beta using loop
 for (i in 1:nrow(indices)) {
   cov_beta[indices[i, 1], indices[i, 2]] <- cov_beta[indices[i, 2], indices[i, 1]] <- Cov_dbeta[i]
 }
 
 cov_VB_beta <- solve(-Hessian_manu - cov_beta)
 time2 <- Sys.time()
 
 results$VB.mean <- output.VB$mu_beta
 if (sum(is.na(sqrt(diag(cov_VB_beta)))) > 0) {
   print(paste0('Seed ',seed, ": NaNs produced"))}
 results$VB.std <- sqrt(diag(output.VB$var_beta))
 results$VB.lower <- output.VB$mu_beta - 1.96*results$VB.std
 results$VB.upper <- output.VB$mu_beta + 1.96*results$VB.std
 results$Coverage.VB <- results$VB.lower <= param & results$VB.upper >= param
 results$VB.std_cor <- sqrt(abs(diag(cov_VB_beta)))
 results$VB.lower_cor <- output.VB$mu_beta - 1.96*results$VB.std_cor
 results$VB.upper_cor <- output.VB$mu_beta + 1.96*results$VB.std_cor
 results$Coverage.VB_cor <- results$VB.lower_cor <= param & results$VB.upper_cor >= param
 results$time.VB <- as.numeric(difftime(time2, time1, units = "secs"))
 results$seed <- seed
 
 return(results)
}

output <- mclapply(((nm_set-1)*n_sim+1):(nm_set*n_sim), 
                   simu_func, mc.cores = nm_cores)
time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
saveRDS(output, file = paste0('VB_MPT_known_', time,'.rds'))


