library("parallel")
library("groupTesting")

Rcpp::sourceCpp('VB_MPT.cpp')

# Compute the log-determinant of a matrix
ldet <- function(X) {
  if(!is.matrix(X)) return(log(X))
  determinant(X,logarithm = TRUE)$modulus
}

########## VB - function: Beta ###############
VB_beta <- function(input_x, input_y, input_omega, input_prior){
  var_beta_inv <- crossprod(input_x*input_omega, input_x) + input_prior$Sigma_inv
  var_beta <- solve(var_beta_inv)
  m_ome <- crossprod(input_x, input_y-1/2)
  mu_beta <- var_beta %*% (m_ome + input_prior$mu_Sigma)
  return(cbind(mu_beta,var_beta))
}

cal_prob0 <- function(yvals, x_abc_input_beta, Se1, Sp1){
  pool_abc <- sum(yvals)
  p1 <- exp((pool_abc>0)*log(1-Se1) + (1-(pool_abc>0))*log(Sp1) 
            + sum(yvals*(x_abc_input_beta))) 
  return(p1)
}

cal_prob1 <- function(yvals, x_abc_input_beta, input_Se, input_Sp){
  pool_abc <- sum(yvals)
  Se1 <- input_Se # pool sensitivity
  Sp1 <- input_Sp # pool specificity
  
  p1 <- exp((pool_abc>0)*log(Se1) + (1-(pool_abc>0))*log(1-Sp1) + sum(yvals*(x_abc_input_beta)))
  
  return(p1)
}

cal_density_Y0 <- function(Eyvals, x_abc, input_beta, 
                           Se1, Sp1){
  pool_abc <- 1- prod(1-Eyvals)
  p1 <- (pool_abc)*log(1-Se1) + (1-pool_abc)*log(Sp1) 
  + sum(Eyvals*(x_abc%*%input_beta))
  return(p1)
}

cal_density_Y1 <- function(Eyvals, x_abc, input_beta, Se, Sp){
  Se1 <- Se[1] # pool sensitivity
  Sp1 <- Sp[1] # pool specificity
  
  pool_abc <- 1- prod(1-Eyvals) 
  p1 <- pool_abc*log(Se1) + (1-pool_abc)*log(1-Sp1) + sum(Eyvals*(x_abc%*%input_beta))
  return(p1)
}


## H Manually Calculation - expectation of second derivatives (same as it)
cal_Hessian_it <- function(mean_beta, mean_y, input_x){
  p_est <- plogis(input_x %*% mean_beta)
  pq_est <- p_est*(1-p_est)
  
  p1 <- sum(-pq_est)
  p2 <- sum(-input_x[,2]^2*pq_est)
  p3 <- sum(-input_x[,3]^2*pq_est)
  
  p12 <- sum(-input_x[,2]*pq_est)
  p13 <- sum(-input_x[,3]*pq_est)
  p23 <- sum(-(input_x[,2]*input_x[,3])*pq_est)
  
  Hessian_manu <- matrix(NA, nrow=dim(input_x)[2], ncol=dim(input_x)[2])
  Hessian_manu <- diag(c(p1,p2,p3))
  Hessian_manu[1,2] <- Hessian_manu[2,1] <- p12
  Hessian_manu[1,3] <- Hessian_manu[3,1] <- p13
  Hessian_manu[2,3] <- Hessian_manu[3,2] <- p23
  
  return(Hessian_manu)
}

#Sample size
N <- 300

#Beta
param <- c(-4,2,1)

## Master Pool Testing (MPT) - known Se&Sp
S <- 1 #1-stage, master pool testing
psz <- 10 #pool size
Se_true <- 0.95 #Pool sensitivity
Sp_true <- 0.98 #Pool specificity
#beta prior
tau_sq <- 100
assayID <- 1 #Assays used in master pools
y.mat <- as.matrix(expand.grid(rep(list(0:1),psz)))
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
  
 ## Group testing data generation
 results.list <- list()
 con.ct <- 1
 x <- cbind(1, rnorm(N), rbinom(N,1,0.5))
 colnames(x) <- c("Intercept", "x1", "x2")

 pReg <- plogis(x %*% param)
 # Simulating test responses
 MPT.gtData <- hier.gt.simulation(N, p=pReg, S, psz, Se=Se_true, Sp=Sp_true, assayID)$gtData
 
 J <- N/psz[1]
 pool_index_mtx <- MPT.gtData[c(1:J),c(6:(6+psz-1))]
 data_x <- x
 data_z <- MPT.gtData[c(1:J),1]
 p <- dim(data_x)[2] # p
 size <- psz[1]      # Pool size
 lowerbound <- numeric(maxiter) #initial value of the ELBO in VB
 
 time1 <- Sys.time()
 ######### Part I. VB CAVI algorithm ################
 PSigma <- tau_sq*diag(p) #Prior of Sigma_beta
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
     pool_data <- VB_ytilde0_cpp(x_sub,
                                 mu_beta,
                                 y.mat,
                                 Se_true,
                                 Sp_true)
     
     vq_z[j] <- cal_density_Y0( pool_data, 
                               x_sub, mu_beta,  Se_true, Sp_true)
     y0[pool_index] <- pool_data
   }
   if (data_z[j] == 1){
     pool_data <- VB_ytilde1_cpp(x_sub,
                                 mu_beta,
                                 y.mat,
                                 Se_true,
                                 Sp_true)
     vq_z[j] <- cal_density_Y1(pool_data, 
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
       pool_data <- VB_ytilde0_cpp(x_sub,
                                   mu_beta,
                                   y.mat,
                                   Se_true,
                                   Sp_true)
       vq_z[j] <- cal_density_Y0(pool_data, 
                                 x_sub, mu_beta,  Se_true, Sp_true)
       y0[pool_index] <- pool_data
     }
     if (data_z[j] == 1){
       pool_data <- VB_ytilde1_cpp(x_sub,
                                   mu_beta,
                                   y.mat,
                                   Se_true,
                                   Sp_true)
       vq_z[j] <- cal_density_Y1(pool_data, 
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
     ## H Manually Calculation
     Hessian_manu <- cal_Hessian_it(output.VB$mu_beta, output.VB$est_y, data_x)
     ## Cov of first derivatives_beta calculation
     Cov_sums <- cal_cov_sums_mpt_cpp(output.VB$mu_beta, output.VB$est_y, 
                                   data_x, data_z, Se_true, Sp_true, pool_index_mtx, y.mat)
     
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
     
     results <- NULL
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
     results.list[[con.ct]] <- results
     
     print(con_tol)
     if(con.ct < length(con_tol_vec)){
      con.ct = con.ct + 1
      con_tol <- con_tol_vec[con.ct]
     }else{
       break
       }
   }
  print(t)
 }
 return(results.list)
}

out.mat <- NULL
time.vec <- NULL
for(sims in 1:n_sim){
  output.cmb = simu_func(1)
  for(j in 1:length(output.cmb)){
    output = output.cmb[[j]]
    out.mat <- rbind(out.mat,output[[1]], output[[2]], output[[3]], output[[4]], output[[5]], output[[6]], output[[7]], output[[8]], output[[9]])
    time.vec <- c(time.vec,output[[10]])
  }
}

write.table(out.mat,'betares.txt',row.names = FALSE, col.names = FALSE)
write.table(time.vec,'timeres.txt',row.names = FALSE, col.names = FALSE)

