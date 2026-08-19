library("parallel")
library("groupTesting")
library(Matrix)
set.seed(1)

#setwd('C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Tick Surveillance R01/Statistical Methods/Variational Bayes')

Rcpp::sourceCpp('C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Tick Surveillance R01/Statistical Methods/Variational Bayes/Data Application/VB_MPT/VBRccP_DA.cpp')

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

cal_prob0.ap <- function(index, x_abc, input_beta, Se1, Sp1){
  yvals <- y.mat[index,]
  pool_abc <- sum(yvals)
  p1 <- exp((pool_abc>0)*log(1-Se1) + (1-(pool_abc>0))*log(Sp1) 
            + sum(yvals*(x_abc%*%input_beta))) 
  return(p1)
}

# VB_ytilde0 <- function(input_x_sub, input_mu_beta,
#                        input_Se, input_Sp){
#   x_sub <- input_x_sub
#   Se1 <- input_Se
#   Sp1 <- input_Sp
#   
#   pool.size <- dim(x_sub)[1]
# 
#   ## Calculate p_abc for all 8 possible combinations under a pool size of 3
#   
#   x_sub_input_mu_beta <- x_sub%*%input_mu_beta
#   p_sub <- rep(NA,dim(y.mat)[1])
#   for(config in 1:length(p_sub)){
#     p_sub[config] <- cal_prob0(y.mat[config,], x_sub_input_mu_beta, Se1, Sp1)
#   }
#   #p_sub <- apply(matrix(seq(from = 1,dim(y.mat)[1])),1,cal_prob0.ap,x_abc = x_sub, input_beta = input_mu_beta, Se1 = Se1, Sp1= Sp1)
#   
#   p_denominator <- sum(p_sub)
#   
#   p_cal <- rep(NA,pool.size)
#   for(ind in 1:pool.size){
#     row.index = which(y.mat[,ind] == 1)
#     #row.index <- row.index.mat[ind,]
#     p_cal[ind] = sum(p_sub[row.index])/p_denominator
#   }
#   
#   return(p_cal)
# }


cal_prob1 <- function(yvals, x_abc_input_beta, input_Se, input_Sp){
  pool_abc <- sum(yvals)
  Se1 <- input_Se # pool sensitivity
  Sp1 <- input_Sp # pool specificity
  
  p1 <- exp((pool_abc>0)*log(Se1) + (1-(pool_abc>0))*log(1-Sp1) + sum(yvals*(x_abc_input_beta)))
  
  return(p1)
}

# VB_ytilde1 <- function(input_x_sub, input_mu_beta, 
#                        input_Se, input_Sp){
#   x_sub <- input_x_sub                     
#   Se1 <- input_Se # pool sensitivity
#   Sp1 <- input_Sp # pool specificity
#   
#   pool.size <- dim(x_sub)[1]
#   
#   ## Calculate p_abc for all 8 possible combinations under a pool size of 3
#   
#   x_sub_input_mu_beta <- x_sub%*%input_mu_beta
#   p_sub <- rep(NA,dim(y.mat)[1])
#   for(config in 1:length(p_sub)){
#     p_sub[config] <- cal_prob1(y.mat[config,], x_sub_input_mu_beta, Se1, Sp1)
#   }
#   
#   p_denominator <- sum(p_sub)
#   
#   p_cal <- rep(NA,pool.size)
#   for(ind in 1:pool.size){
#     row.index = which(y.mat[,ind] == 1)
#     #row.index <- row.index.mat[ind,]
#     p_cal[ind] = sum(p_sub[row.index])/p_denominator
#   }
#   
#   return(p_cal)
# }

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


# ## H Manually Calculation - expectation of second derivatives (same as it)
# cal_Hessian_it <- function(mean_beta, mean_y, input_x){
#   p_est <- plogis(input_x %*% mean_beta)
#   pq_est <- p_est*(1-p_est)
#   
#   p1 <- sum(-pq_est)
#   p2 <- sum(-input_x[,2]^2*pq_est)
#   p3 <- sum(-input_x[,3]^2*pq_est)
#   
#   p12 <- sum(-input_x[,2]*pq_est)
#   p13 <- sum(-input_x[,3]*pq_est)
#   p23 <- sum(-(input_x[,2]*input_x[,3])*pq_est)
#   
#   Hessian_manu <- matrix(NA, nrow=dim(input_x)[2], ncol=dim(input_x)[2])
#   Hessian_manu <- diag(c(p1,p2,p3))
#   Hessian_manu[1,2] <- Hessian_manu[2,1] <- p12
#   Hessian_manu[1,3] <- Hessian_manu[3,1] <- p13
#   Hessian_manu[2,3] <- Hessian_manu[3,2] <- p23
#   
#   return(Hessian_manu)
# }

## H Manually Calculation - expectation of second derivatives (same as it)
cal_Hessian_it <- function(mean_beta, mean_y, input_x){
  size = psz[1]
  p_est <- plogis(input_x %*% mean_beta)
  pq_est <- p_est*(1-p_est)
  
  Hessian_manu <- matrix(NA, nrow=dim(input_x)[2], ncol=dim(input_x)[2])
  for(ii in 1:dim(Hessian_manu)[1]){
    for(jj in 1:ii){
      Hessian_manu[ii,jj] = sum(-input_x[,ii]*input_x[,jj]*pq_est)
      Hessian_manu[jj,ii] = sum(-input_x[,ii]*input_x[,jj]*pq_est)
    }
  }
  
  return(Hessian_manu)
}

## Summations of Covariance_y1y2 * multipliers
# cal_cov_sums_mpt <- function(mean_beta, mean_y,
#                              input_x, input_z, input_Se, input_Sp,
#                              input_pool_index){
#   J <- dim(input_pool_index)[1]
#   size <- dim(input_pool_index)[2]
# 
#   # Calculate cov(y_i,y_j) = E(y_i*y_j) - E(y_i)E(y_j)
#   Eyy_cov <- matrix(NA, J, size*(size-1)/2)
# 
#   Eyy_cov_sum11 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta0
#   Eyy_cov_sum22 <- matrix(NA, J, size*(size-1)/2) #beta1 vs. beta1
#   Eyy_cov_sum33 <- matrix(NA, J, size*(size-1)/2) #beta2 vs. beta2
# 
#   Eyy_cov_sum1 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta1
#   Eyy_cov_sum2 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta2
#   Eyy_cov_sum3 <- matrix(NA, J, size*(size-1)/2) #beta1 vs. beta2
# 
#   for (j in 1:J){
#     pool_index <- input_pool_index[j,]
#     x_sub <- input_x[pool_index,]
#     Ey_spool <- mean_y[pool_index]
#     x_sub_mean_beta <- x_sub%*%mean_beta
#     if (input_z[j] == 0){
#       p_sub <- apply(y.mat, 1,
#                      function(row) cal_prob0(row,
#                                              x_sub_mean_beta, input_Se, input_Sp))
#       p_denominator <- sum(p_sub)
#     }
#     if (input_z[j] == 1){
#       p_sub <- apply(y.mat, 1,
#                      function(row) cal_prob1(row,
#                                              x_sub_mean_beta, input_Se, input_Sp))
#       p_denominator <- sum(p_sub)
#     }
#     # Eyy.12 <- sum(p_sub[4], p_sub[8])/p_denominator #E(y1*y2)
#     # Eyy.13 <- sum(p_sub[6], p_sub[8])/p_denominator #E(y1*y3)
#     # Eyy.23 <- sum(p_sub[7], p_sub[8])/p_denominator #E(y2*y3)
# 
#     Eyy <- rep(NA,size*(size-1)/2)
#     ct = 1
#     for(ii in 1:(size-1)){
#       for(jj in (ii+1):size){
#         index = y.mat[,ii] == 1 & y.mat[,jj] == 1
#         Eyy[ct] = sum(p_sub[index])/p_denominator
#         ct = ct + 1
#       }
#     }
# 
#     # Eyy_cov[j,1] <- Eyy.12 - Ey_spool[1]*Ey_spool[2] #E(y1*y2)-E(y1)E(y2)
#     # Eyy_cov[j,2] <- Eyy.13 - Ey_spool[1]*Ey_spool[3] #E(y1*y3)-E(y1)E(y3)
#     # Eyy_cov[j,3] <- Eyy.23 - Ey_spool[2]*Ey_spool[3] #E(y2*y3)-E(y2)E(y3)
# 
#     ct = 1
#     for(ii in 1:(size-1)){
#       for(jj in (ii+1):size){
#         Eyy_cov[j,ct] = Eyy[ct] - Ey_spool[ii]*Ey_spool[jj]
#         ct = ct + 1
#       }
#     }
# 
#     # x_sub_multi <- matrix(c(x_sub[1,]*x_sub[2,] + x_sub[2,]*x_sub[1,],
#     #                         x_sub[1,]*x_sub[3,] + x_sub[3,]*x_sub[1,],
#     #                         x_sub[2,]*x_sub[3,] + x_sub[3,]*x_sub[2,]), nrow = size)
# 
#     x_sub_multi <- matrix(NA,(size*(size-1)/2),dim(x_sub)[2])
#     ct = 1
#     for(ii in 1:(size-1)){
#       for(jj in (ii+1):size){
#         x_sub_multi[ct,] = x_sub[ii,]*x_sub[jj,] + x_sub[jj,]*x_sub[ii,]
#         ct = ct + 1
#       }
#     }
#     x_sub_multi <- t(x_sub_multi)
# 
#     Eyy_cov_sum11[j,] <- x_sub_multi[1,] * Eyy_cov[j,]
#     Eyy_cov_sum22[j,] <- x_sub_multi[2,] * Eyy_cov[j,]
#     Eyy_cov_sum33[j,] <- x_sub_multi[3,] * Eyy_cov[j,]
# 
#     # Eyy_cov_sum1[j,1] <- (x_sub[1,2]+x_sub[2,2])*Eyy_cov[j,1]
#     # Eyy_cov_sum1[j,2] <- (x_sub[1,2]+x_sub[3,2])*Eyy_cov[j,2]
#     # Eyy_cov_sum1[j,3] <- (x_sub[2,2]+x_sub[3,2])*Eyy_cov[j,3]
#     #
#     # Eyy_cov_sum2[j,1] <- (x_sub[1,3]+x_sub[2,3])*Eyy_cov[j,1]
#     # Eyy_cov_sum2[j,2] <- (x_sub[1,3]+x_sub[3,3])*Eyy_cov[j,2]
#     # Eyy_cov_sum2[j,3] <- (x_sub[2,3]+x_sub[3,3])*Eyy_cov[j,3]
#     #
#     # Eyy_cov_sum3[j,1] <- (x_sub[1,2]*x_sub[2,3]+x_sub[2,2]*x_sub[1,3])*Eyy_cov[j,1]
#     # Eyy_cov_sum3[j,2] <- (x_sub[1,2]*x_sub[3,3]+x_sub[3,2]*x_sub[1,3])*Eyy_cov[j,2]
#     # Eyy_cov_sum3[j,3] <- (x_sub[2,3]*x_sub[3,3]+x_sub[3,2]*x_sub[2,3])*Eyy_cov[j,3]
# 
#     ct = 1
#     for(ii in 1:(size-1)){
#       for(jj in (ii+1):size){
#         Eyy_cov_sum1[j,ct] <- (x_sub[ii,2]+x_sub[jj,2])*Eyy_cov[j,ct]
#         Eyy_cov_sum2[j,ct] <- (x_sub[ii,3]+x_sub[jj,3])*Eyy_cov[j,ct]
#         Eyy_cov_sum3[j,ct] <- (x_sub[ii,2]*x_sub[jj,3]+x_sub[jj,2]*x_sub[ii,3])*Eyy_cov[j,ct]
#         ct = ct + 1
#       }
#     }
# 
# 
#   }
#   return(c(sum(Eyy_cov_sum11), sum(Eyy_cov_sum22), sum(Eyy_cov_sum33),
#            sum(Eyy_cov_sum1), sum(Eyy_cov_sum2), sum(Eyy_cov_sum3)))
# 
# }

# ## Summations of Covariance_y1y2 * multipliers 
# cov_Y_mpt <- function(mean_beta, mean_y, 
#                              input_x, input_z, input_Se, input_Sp,
#                              input_pool_index){
#   J <- dim(input_pool_index)[1]
#   size <- dim(input_pool_index)[2]
#   
#   # Calculate cov(y_i,y_j) = E(y_i*y_j) - E(y_i)E(y_j)
#   Eyy_cov <- matrix(NA, J, size*(size-1)/2)
#   
#   Eyy_cov_sum11 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta0
#   Eyy_cov_sum22 <- matrix(NA, J, size*(size-1)/2) #beta1 vs. beta1
#   Eyy_cov_sum33 <- matrix(NA, J, size*(size-1)/2) #beta2 vs. beta2
#   
#   Eyy_cov_sum1 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta1
#   Eyy_cov_sum2 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta2
#   Eyy_cov_sum3 <- matrix(NA, J, size*(size-1)/2) #beta1 vs. beta2
#   
#   for (j in 1:J){
#     pool_index <- input_pool_index[j,]
#     x_sub <- input_x[pool_index,]
#     Ey_spool <- mean_y[pool_index]
#     x_sub_mean_beta <- x_sub%*%mean_beta
#     if (input_z[j] == 0){
#       p_sub <- apply(y.mat, 1, 
#                      function(row) cal_prob0(row,
#                                              x_sub_mean_beta, input_Se, input_Sp))
#       p_denominator <- sum(p_sub)
#     }
#     if (input_z[j] == 1){
#       p_sub <- apply(y.mat, 1, 
#                      function(row) cal_prob1(row,
#                                              x_sub_mean_beta, input_Se, input_Sp))
#       p_denominator <- sum(p_sub)
#     }
#     
#     Eyy <- rep(NA,size*(size-1)/2)
#     ct = 1
#     for(ii in 1:(size-1)){
#       for(jj in (ii+1):size){
#         index = y.mat[,ii] == 1 & y.mat[,jj] == 1
#         Eyy[ct] = sum(p_sub[index])/p_denominator
#         ct = ct + 1
#       }
#     }
#     
#     ct = 1
#     for(ii in 1:(size-1)){
#       for(jj in (ii+1):size){
#         Eyy_cov[j,ct] = Eyy[ct] - Ey_spool[ii]*Ey_spool[jj]
#         ct = ct + 1
#       }
#     }
#   }
#     return(Eyy_cov)
# }



##### Import Data ######
path = 'C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Tick Surveillance R01/Statistical Methods/Variational Bayes/Data Application/CleanedData/CleanedData/'
psz = 3
data_z <- readRDS(paste0(path,"Data_z_pz3.rds"))
J = length(data_z)
N = J*psz
data_x <- readRDS(paste0(path,"Data_x.rds"))
data_x <- data_x[1:N,]

#retest_y <- readRDS(paste0(path,"Retest_y_pz3.rds"))
#retest_y_sub <- readRDS(paste0(path,"Retest_y_sub_pz3.rds"))
#retest_ind <-readRDS(paste0(path,"Retest_ind_pz3.rds"))
pool_index_mtx <- readRDS(paste0(path,"pool_index_mtx_pz3.rds"))

## Master Pool Testing (MPT) - known Se&Sp
Se_true <- 0.9734 #Pool sensitivity
Sp_true <- 0.9999 #Pool specificity
#beta prior
tau_sq <- 100
assayID <- 1 #Assays used in master pools
y.mat <- as.matrix(expand.grid(rep(list(0:1),psz)))
row.index.mat <- matrix(NA,dim(y.mat)[2],dim(y.mat)[1])
for(r in 1:dim(y.mat)[2]){
  row.index.mat[r,] = y.mat[,r] == 1
}
p = dim(data_x)[2]

# Define the indices for combinations
indices <- t(combn(c(1:p),2))

con_tol_vec <- c(1,0.1,0.01,10^-3,10^-4,10^-5,10^-6,10^-7,10^-8) #converge criteria
time.vec <- rep(NA,length(con_tol_vec)+1)
con_tol <- con_tol_vec[1]
maxiter <- 10000 #maximum iterations in VB

results.list <- list()
con.ct <- 1
lowerbound <- numeric(maxiter) #initial value of the ELBO in VB
time.vec[1] <- Sys.time()
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
pool_data <- rep(NA, psz)
for (j in 1:J){
  pool_index <- pool_index_mtx[j,]
  x_sub <- data_x[pool_index,]
  if (data_z[j] == 0){
    pool_data <- VB_ytilde0_cpp(x_sub, mu_beta, y.mat, Se_true, Sp_true)
    vq_z[j] <- cal_density_Y0( pool_data, 
                               x_sub, mu_beta,  Se_true, Sp_true)
    y0[pool_index] <- pool_data
  }
  if (data_z[j] == 1){
    pool_data <- VB_ytilde1_cpp(x_sub, mu_beta,y.mat, Se_true, Sp_true)
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
     pool_data <- VB_ytilde0_cpp(x_sub, mu_beta,y.mat, Se_true, Sp_true)
     vq_z[j] <- cal_density_Y0(pool_data, 
                               x_sub, mu_beta,  Se_true, Sp_true)
     y0[pool_index] <- pool_data
   }
   if (data_z[j] == 1){
     pool_data <- VB_ytilde1_cpp(x_sub, mu_beta,y.mat, Se_true, Sp_true)
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
   #Cov_sums <- cal_cov_sums_mpt(output.VB$mu_beta, output.VB$est_y, 
   #                             data_x, data_z, Se_true, Sp_true, pool_index_mtx)
     
   Ey_1mins <- output.VB$est_y * (1-output.VB$est_y)
   #Var_dbeta <- colSums(data_x^2*Ey_1mins) + Cov_sums[1:p] #Variance of first derivatives: beta0, beta1, beta2
   cov.Y.vals <- cov_Y_mpt_cpp(output.VB$mu_beta, output.VB$est_y, 
                           data_x, data_z, Se_true, Sp_true, pool_index_mtx,y.mat)
   # cov.Y.mat <- sparseMatrix(i = 1, j = 1,x = 0, dims = c(N,N))
   # for(j in 1:J){
   #   sub.mat <- matrix(0,psz,psz)
   #   diag(sub.mat) <- Ey_1mins[((j-1)*psz + 1):(j*psz)]
   #   ct = 1
   #   for(ii in 2:psz){
   #     for(jj in (1:(ii-1))){
   #       sub.mat[ii,jj] = cov.Y.vals[j,ct]
   #       sub.mat[jj,ii] = cov.Y.vals[j,ct]
   #       ct = ct + 1
   #     }
   #   }
   #   cov.Y.mat[((j-1)*psz + 1):(j*psz), ((j-1)*psz + 1):(j*psz)] <- sub.mat
   # }
   cov_beta <- cov_beta_cpp(
     data_x,
     cov.Y.vals,
     Ey_1mins,
     psz
   )
   # # Covariance of first derivatives: beta0-1, beta0-2, beta1-2
   #   multipliers <- apply(indices, 1, 
   #                        function(vec) data_x[,vec[1]] * data_x[, vec[2]])
   #   Cov_dbeta <- colSums(multipliers*Ey_1mins) + Cov_sums[c(p+1):(2*p)]
   #   
   #   cov_beta <- diag(Var_dbeta)
   #   # Assign values from cov_xy to cov_beta using loop
   #   for (i in 1:nrow(indices)) {
   #     cov_beta[indices[i, 1], indices[i, 2]] <- cov_beta[indices[i, 2], indices[i, 1]] <- Cov_dbeta[i]
   #   }
   #   
     cov_VB_beta <- solve(-Hessian_manu - cov_beta)
     time.vec[con.ct+1] <- Sys.time()
     
     results <- NULL
     results$VB.mean <- output.VB$mu_beta
     if (sum(is.na(sqrt(diag(cov_VB_beta)))) > 0) {
       print(paste0('Seed ',seed, ": NaNs produced"))}
     results$VB.std <- sqrt(diag(output.VB$var_beta))
     results$VB.lower <- output.VB$mu_beta - 1.96*results$VB.std
     results$VB.upper <- output.VB$mu_beta + 1.96*results$VB.std
     results$VB.std_cor <- sqrt(abs(diag(cov_VB_beta)))
     results$VB.lower_cor <- output.VB$mu_beta - 1.96*results$VB.std_cor
     results$VB.upper_cor <- output.VB$mu_beta + 1.96*results$VB.std_cor
     #results$seed <- seed
     results.list[[con.ct]] <- results
     
     print(con_tol)
     if(con.ct < length(con_tol_vec)){
      con.ct = con.ct + 1
      con_tol <- con_tol_vec[con.ct]
     }else{
       break
       }
   }
  #print(t)
 }

# for(i in 2:length(con_tol_vec)){
#   results1 = results.list[[i-1]]
#   results2 = results.list[[i]]
#   delta <- results1$VB.mean - results2$VB.mean
#   if(max(abs(delta))<0.001){
#     conv.iter = i
#     time = time.vec[i]-time.vec[1]
#     break
#   }
# }
# 
# results.list[[conv.iter]]
# time

conv.iter = 4
con_tol_vec[conv.iter]
results.list[[conv.iter]]
results.list[[conv.iter]]
out.df = data.frame(results.list[[conv.iter]][1],results.list[[conv.iter]][2],results.list[[conv.iter]][3],results.list[[conv.iter]][4],results.list[[conv.iter]][5],results.list[[conv.iter]][6],results.list[[conv.iter]][7] )
out.df = round(out.df,3)
amp.col = rep("&", dim(out.df)[1])
print.df = data.frame(out.df[,1], amp.col,out.df[,2], amp.col,rep("(", dim(out.df)[1]), out.df[,3], rep(",", dim(out.df)[1]),out.df[,4],rep(")", dim(out.df)[1]), out.df[,5], amp.col,rep("(", dim(out.df)[1]), out.df[,6], rep(",", dim(out.df)[1]), out.df[,7], rep(")", dim(out.df)[1]),amp.col )
names(print.df) = seq(1,dim(print.df)[2])
print.df

time.vec[conv.iter + 1] - time.vec[1]
#save.image('C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Tick Surveillance R01/Statistical Methods/Variational Bayes/Data Application/VB_MPT/VB_MPT_DA.RData')
