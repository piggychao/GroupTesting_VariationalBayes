library("parallel")
library("groupTesting")
set.seed(1)

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

######### VB - function: Y_tilde ################
### when data_z=0
cal_prob0 <- function(a, x_a, input_beta, Se1, Sp1){
  
  p1 <- exp((a>0)*log(1-Se1) + (1-(a>0))*log(Sp1) + 
              a*(x_a %*% input_beta)) 
  return(p1)
}

VB_ytilde0 <- function(input_x_sub, input_mu_beta,
                       input_Se, input_Sp){
  x_a <- input_x_sub
  
  ## Calculate p_a for all 2 possible combinations under a pool size of 1
  p_sub1 <- cal_prob0(0, x_a, input_mu_beta, input_Se, input_Sp)
  p_sub2 <- cal_prob0(1, x_a, input_mu_beta, input_Se, input_Sp)
  
  p_denominator <- sum(p_sub1, p_sub2)
  
  p_cal1 <- p_sub2/p_denominator
  
  return(p_cal1)
}

### when data_z=1 
cal_prob1 <- function(a, x_a, input_beta, input_Se, input_Sp){
  Se1 <- input_Se # pool sensitivity
  Sp1 <- input_Sp # pool specificity
  
  p1 <- exp((a>0)*log(Se1) + (1-(a>0))*log(1-Sp1) + 
              a*(x_a %*% input_beta))
  return(p1)
}

VB_ytilde1 <- function(input_x_sub, input_mu_beta, 
                       input_Se, input_Sp){
  x_a <- input_x_sub
  
  ## Calculate p_a for all 2 combinations under a pool size of 1
  p_sub1 <- cal_prob1(0, x_a, input_mu_beta, input_Se, input_Sp)
  p_sub2 <- cal_prob1(1, x_a, input_mu_beta, input_Se, input_Sp)
  
  p_denominator <- sum(p_sub1, p_sub2)
  
  p_cal1 <- p_sub2/p_denominator
  
  return(p_cal1)
}

############ Calculate Y_tilde Joint Density Function #####
cal_density_Y0 <- function(Ey, x_a, input_beta, Se1, Sp1){
  p1 <- Ey*log(1-Se1) + (1-Ey)*log(Sp1) + 
    Ey*(x_a %*% input_beta) 
  return(p1)
}

cal_density_Y1 <- function(Ey, x_a, input_beta, Se1, Sp1){
  p1 <- Ey*log(Se1) + (1-Ey)*log(1-Sp1) + 
    Ey*(x_a %*% input_beta) 
  return(p1)
}

# ## H Manually Calculation - expectation of second derivatives
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

## Covariance of first derivatives
cal_covar_it <- function(mean_y, input_x){
  # Define the indices for combinations
  indices <- t(combn(c(1:dim(input_x)[2]),2))
  
  Ey1min <- mean_y * (1-mean_y)
  
  var_y <- colSums(input_x^2 * Ey1min)
  cov_d0d1 <- sum(input_x[,1]*input_x[,2]*Ey1min) #d1: beta0 vs beta1
  cov_d0d2 <- sum(input_x[,1]*input_x[,3]*Ey1min) #d1: beta0 vs beta2
  cov_d1d2 <- sum(input_x[,2]*input_x[,3]*Ey1min) #d1: beta1 vs beta2
  cov_xy <- c(cov_d0d1, cov_d0d2, cov_d1d2)
  
  cov_beta <- diag(var_y)
  # Assign values from cov_xy to cov_beta using loop
  for (i in 1:nrow(indices)) {
    cov_beta[indices[i, 1], indices[i, 2]] <- cov_beta[indices[i, 2], indices[i, 1]] <- cov_xy[i]
  }
  
  return(cov_beta)
}  

##### Import Data ######
path = 'C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Tick Surveillance R01/Statistical Methods/Variational Bayes/Data Application/CleanedData/CleanedData/'
psz = 1
data_z <- readRDS(paste0(path,"Data_z_pz1.rds"))
J = length(data_z)
N = J*psz
data_x <- readRDS(paste0(path,"Data_x.rds"))
data_x <- data_x[1:N,]


## Individual Pool Testing (IT) - known Se&Sp
S <- 1 #1-stage, individual pool testing 
psz <- 1 #pool size
Se_true <- 0.9734 #Pool sensitivity
Sp_true <- 0.9999 #Pool specificity
assayID <- 1 #Assays used in master pools
tau_sq <- 1

con_tol_vec <- c(1,0.1,0.01,10^-3,10^-4,10^-5,10^-6,10^-7,10^-8) #converge criteria
con_tol <- con_tol_vec[1]
maxiter <- 10000 #maximum iterations in VB

results.list <- list()
con.ct <- 1
p <- dim(data_x)[2] # p
size <- psz[1]      # Pool size
lowerbound <- numeric(maxiter) #initial value of the ELBO in VB

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
 x_sub <- data_x[j,]
 if (data_z[j] == 0){
   pool_data <- VB_ytilde0(x_sub, mu_beta, Se_true, Sp_true)
   vq_z[j] <- cal_density_Y0(pool_data, x_sub, mu_beta,  Se_true, Sp_true)
   y0[j] <- pool_data
 }
 if (data_z[j] == 1){
   pool_data <- VB_ytilde1(x_sub, mu_beta, Se_true, Sp_true)
   vq_z[j] <- cal_density_Y1(pool_data,x_sub, mu_beta, Se_true, Sp_true)
   y0[j] <- pool_data
 }
}
 
lowerbound[1] <- sum(vq_z) + 
 0.5*p + 0.5*ldet(var_beta) + 0.5*Pdet - 0.5*t(mu_beta - Pmu)%*%PSigma_inv%*%(mu_beta - Pmu) + 
 sum((y0-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(PSigma_inv %*% var_beta))

# Iteration t=2: maxiter
output.VB <- list()
time.vec <- as.numeric(Sys.time())
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
   x_sub <- data_x[j,]
   if (data_z[j] == 0){
     pool_data <- VB_ytilde0(x_sub, mu_beta, Se_true, Sp_true)
     vq_z[j] <- cal_density_Y0(pool_data, x_sub, mu_beta,  Se_true, Sp_true)
     y0[j] <- pool_data
   }
   if (data_z[j] == 1){
     pool_data <- VB_ytilde1(x_sub, mu_beta, Se_true, Sp_true)
     vq_z[j] <- cal_density_Y1(pool_data, x_sub, mu_beta, Se_true, Sp_true)
     y0[j] <- pool_data
   }
  }
  
  lowerbound[t]  <- sum(vq_z) + 
   0.5*p + 0.5*ldet(var_beta) + 0.5*Pdet - 0.5*t(mu_beta - Pmu)%*%PSigma_inv%*%(mu_beta - Pmu) + 
   sum((y0-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(PSigma_inv %*% var_beta))
   
  if(abs(lowerbound[t] - lowerbound[t-1]) < con_tol) {
   output.VB <- list(mu_beta=mu_beta, 
                     var_beta=var_beta, 
                     est_y=y0)
     Hessian_p1 <- cal_Hessian_it(output.VB$mu_beta, output.VB$est_y, data_x) #VB_mean_beta, 
     #Hessian_p2 <- cal_covar_it(output.VB$est_y, data_x)
     Hessian_p2 = t(data_x)%*%((output.VB$est_y*(1-output.VB$est_y))*data_x)
     
     cov_VB_beta <- solve(- Hessian_p1 - Hessian_p2)
     
     results <- NULL
     results$VB.mean <- output.VB$mu_beta
     # if (sum(is.na(sqrt(diag(cov_VB_beta)))) > 0) {
     #   print(paste0('Seed ',seed, ": NaNs produced"))}
     results$VB.std <- sqrt(diag(output.VB$var_beta))
     results$VB.lower <- output.VB$mu_beta - 1.96*results$VB.std
     results$VB.upper <- output.VB$mu_beta + 1.96*results$VB.std
     results$VB.std_cor <- sqrt(abs(diag(cov_VB_beta)))
     results$VB.lower_cor <- output.VB$mu_beta - 1.96*results$VB.std_cor
     results$VB.upper_cor <- output.VB$mu_beta + 1.96*results$VB.std_cor
     results.list[[con.ct]] <- results
     time.vec <- c(time.vec, as.numeric(Sys.time()))
     print(con_tol)
     if(con.ct < length(con_tol_vec)){
       con.ct = con.ct + 1
       con_tol <- con_tol_vec[con.ct]
     }else{
       break
     }
   }
   
   #print(paste0(t, ": ",lowerbound[t], " , ", abs(lowerbound[t] - lowerbound[t-1])))
 }

###Convergence check
for(i in 2:length(con_tol_vec)){
  results1 = results.list[[i-1]]
  results2 = results.list[[i]]
  delta <- results1$VB.mean - results2$VB.mean
  if(max(abs(delta))<0.001){
    conv.iter = i
    time = time.vec[i]-time.vec[1]
    break
  }
}

conv.iter = 4
con_tol_vec[conv.iter]
time.vec[conv.iter + 1] - time.vec[1]
results.list[[conv.iter]]
out.df = data.frame(results.list[[conv.iter]][1],results.list[[conv.iter]][2],results.list[[conv.iter]][3],results.list[[conv.iter]][4],results.list[[conv.iter]][5],results.list[[conv.iter]][6],results.list[[conv.iter]][7] )
out.df = round(out.df,3)
amp.col = rep("&", dim(out.df)[1])
print.df = data.frame(out.df[,1], amp.col,out.df[,2], amp.col,rep("(", dim(out.df)[1]), out.df[,3], rep(",", dim(out.df)[1]),out.df[,4],rep(")", dim(out.df)[1]), amp.col, out.df[,5], amp.col,rep("(", dim(out.df)[1]), out.df[,6], rep(",", dim(out.df)[1]), out.df[,7], rep(")", dim(out.df)[1]),amp.col )
names(print.df) = seq(1,dim(print.df)[2])
print.df
save.image('C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Tick Surveillance R01/Statistical Methods/Variational Bayes/Data Application/VB_IT/VB_IT_DA.RData')
