library("parallel")
library("groupTesting")

Rcpp::sourceCpp('DA_VB_DT.cpp')

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
  
  Hessian_manu <- matrix(0, nrow=(p+4), ncol=(p+4))
  for(ii in 1:p){
    for(jj in 1:ii){
      Hessian_manu[ii,jj] = sum(-input_x[,ii]*input_x[,jj]*pq_est)
      Hessian_manu[jj,ii] = sum(-input_x[,ii]*input_x[,jj]*pq_est)
    }
  }
  
  p4 <- x_x1mins(input_Se[1])^2 * sum(-Ez*(input_Se[1]^2+input_z-2*input_z*input_Se[1])) #d2_Se_p
  p5 <- x_x1mins(input_Se[2])^2 * sum(-Ey_sub*(input_Se[2]^2+retest_y_sub-2*retest_y_sub*input_Se[2])) #d2_Se_i
  p6 <- x_x1mins(input_Sp[1])^2 * sum(-Ez_1mins*(input_Sp[1]^2+(1-input_z)-2*(1-input_z)*input_Sp[1])) #d2_Sp_p
  p7 <- x_x1mins(input_Sp[2])^2 * sum(-Ey_sub_1mins*(input_Sp[2]^2+(1-retest_y_sub)-2*(1-retest_y_sub)*input_Sp[2])) #d2_Sp_i
  
  Hessian_manu[(p+1),(p+1)] = p4
  Hessian_manu[(p+2),(p+2)] = p5
  Hessian_manu[(p+3),(p+3)] = p6
  Hessian_manu[(p+4),(p+4)] = p7
  
  return(Hessian_manu)
}

cal_cov_sums_unknw <- function(mean_beta, mean_y,
                               input_x, input_z, input_retest_y,
                               input_Se, input_Sp,
                               input_Se_alp_beta, input_Sp_alp_beta,
                               input_pool_index,mean_z){
  J <- length(input_z)
  size <- dim(input_pool_index)[2]
  p <- dim(input_x)[2]
  
  Ez.mtx <- matrix(1-mean_y, nrow = J, byrow=T)
  Ez_1mins <- 1-mean_z
  Ez <- 1-Ez_1mins # E(Z_tilde=1)
  
  Eyy <- rep(NA,size*(size-1)/2)
  
  # Calculate cov(y_i,y_j) = E(y_i*y_j) - E(y_i)E(y_j)
  Eyy_cov <- matrix(NA, J, size*(size-1)/2)
  
  Eyy_cov_sum11 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta0
  Eyy_cov_sum22 <- matrix(NA, J, size*(size-1)/2) #beta1 vs. beta1
  Eyy_cov_sum33 <- matrix(NA, J, size*(size-1)/2) #beta2 vs. beta2
  
  Eyy_cov_sum44 <- rep(0, J) #Var(Se_i)
  Eyy_cov_sum55 <- rep(0, J) #Var(Sp_i)
  
  Eyy_cov_sum1 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta1
  Eyy_cov_sum2 <- matrix(NA, J, size*(size-1)/2) #beta0 vs. beta2
  Eyy_cov_sum3 <- matrix(NA, J, size*(size-1)/2) #beta1 vs. beta2
  
  ## Cov(Se_p, par:beta0, beta1, beta2)
  cov_Se1par_sum <- matrix(NA, nrow=J, ncol=p)
  ## Cov(Sp_p, par:beta0, beta1, beta2
  cov_Sp1par_sum <- matrix(NA, nrow=J, ncol=p)
  ## Cov(Se_i, par:beta0, beta1, beta2)
  cov_Se2par_sum <- matrix(0,J,p)
  ## Cov(Sp_i, par:beta0, beta1, beta2)
  cov_Sp2par_sum <- matrix(0,J,p)

  cov_Sp2Se2_sum <- rep(0, J) #Cov(Se_i vs. Sp_i)
  
  loop.out <- cal_cov_loop_cpp(mean_beta,
    mean_y,
    input_x,
    input_z,
    input_retest_y,
    input_Se,
    input_Sp,
    input_Se_alp_beta,
    input_Sp_alp_beta,
    input_pool_index,
    mean_z,
    y.mat)
  
  Eyy_cov = loop.out$Eyy_cov
  Eyy_cov_sum44 = loop.out$Eyy_cov_sum44
  Eyy_cov_sum55 = loop.out$Eyy_cov_sum55
  cov_Se1par_sum = loop.out$cov_Se1par_sum
  cov_Sp1par_sum = loop.out$cov_Sp1par_sum
  cov_Se2par_sum = loop.out$cov_Se2par_sum
  cov_Sp2par_sum = loop.out$cov_Sp2par_sum
  cov_Sp2Se2_sum = loop.out$cov_Sp2Se2_sum
  return(list(cov_sums = c(sum(Eyy_cov_sum44, na.rm = T), sum(Eyy_cov_sum55, na.rm = T),
           colSums(cov_Se1par_sum), colSums(cov_Sp1par_sum),
           colSums(cov_Se2par_sum,na.rm = T),
           colSums(cov_Sp2par_sum,na.rm = T),
            sum(cov_Sp2Se2_sum, na.rm = T)),
           Eyy_cov = Eyy_cov))
}

psz = 3
data_z <- readRDS("Data_z_pz3.rds")
J = length(data_z)
N = J*psz
data_x <- readRDS("Data_x.rds")
data_x <- data_x[1:N,]
p <- dim(data_x)[2]

retest_y <- readRDS("Retest_y_pz3.rds")
retest_y_sub <- readRDS("Retest_y_sub_pz3.rds")
retest_ind <-readRDS("Retest_ind_pz3.rds")
pool_index_mtx <- readRDS("pool_index_mtx_pz3.rds")

######### Part I. Generate Data ###############
######   (Use 'groupTesting' package)   ######

## Dorfman Pool Testing (DT) - unknown Se&Sp
S <- 2 #2-stage hierarchical testing
psz <- c(psz,1) #Pool sizes used in stages 1-2
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
indices <- t(combn(c(1:p),2))

con_tol_vec <- c(1,0.1,0.01,10^-3)#,10^-4,10^-5,10^-6,10^-7,10^-8) #converge criteria
con_tol <- con_tol_vec[1]
maxiter <- 10000 #maximum iterations in VB

results.list <- list()
con.ct <- 1

############## Part II. VB Data Preparation ############
size <- psz[1] # Pool size

############ CAVI algorithm ################
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
time1 <- Sys.time()
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
  
  print(t)
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

    cov.out =  cal_cov_sums_unknw(results$mu_beta, results$est_y, 
                                   data_x, data_z, retest_y, 
                                   results$Se, results$Sp, 
                                   results$Se_alp_beta, results$Sp_alp_beta, 
                                   pool_index_mtx, results$est_z)
    cov_sums = cov.out$cov_sums
    cov.Y.vals = cov.out$Eyy_cov
    
    Ey_1mins <- results$est_y * (1-results$est_y)
    cov_beta_sub <- cov_beta_cpp(
      data_x,
      cov.Y.vals,
      Ey_1mins,
      psz[1]
    )

      var_Se1 <- (x_x1mins(results$Se[1])^2) * sum((data_z-results$Se[1])^2 * Ez * Ez_1mins)
      var_Sp1 <- (x_x1mins(results$Sp[1])^2) * sum((1-data_z-results$Sp[1])^2 * Ez * Ez_1mins)
      var_Se2 <- (x_x1mins(results$Se[2])^2) * (sum((retest_y_sub-results$Se[2])^2*Ey_sub*Ey_sub_1mins) +
                                                  cov_sums[1])
      var_Sp2 <- (x_x1mins(results$Sp[2])^2) * (sum((1-retest_y_sub-results$Sp[2])^2*Ey_sub*Ey_sub_1mins) +
                                                  cov_sums[2])
      
      cov_par <- matrix(0,(p+4),(p+4))
      for(pp in 1:p){
        cov_par[pp,1:p] <- cov_beta_sub[pp,1:p]
      }
      cov_par[(p+1),(p+1)] <- var_Se1
      cov_par[(p+2),(p+2)] <- var_Se2
      cov_par[(p+3),(p+3)] <- var_Sp1
      cov_par[(p+4),(p+4)] <- var_Sp2
      
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
      
      cov_Se1_par <- x_x1mins(results$Se[1]) * cov_sums[3:(3+p-1)] #d_Se_p vs Beta's
      cov_Sp1_par <- x_x1mins(results$Sp[1]) * cov_sums[(3+p):(3+2*p-1)] #d_Sp_p vs Beta's
      cov_Se2_par <- x_x1mins(results$Se[2]) * cov_sums[(3+2*p):(3+3*p-1)] #d_Se_i vs Beta's
      cov_Sp2_par <- x_x1mins(results$Sp[2]) * cov_sums[(3+3*p):(3+4*p-1)] #d_Sp_i vs Beta's
      
      cov_Se2_Sp2 <- x_x1mins(results$Se[2]) * x_x1mins(results$Sp[2]) *  cov_sums[(3+4*p)] #d_Se_i vs d_Sp_i
      
      # Assign cov_SeSp_par manually
      cov_par[(p+1),c(1:p)] <- cov_par[c(1:p),(p+1)] <- cov_Se1_par
      cov_par[(p+2),c(1:p)] <- cov_par[c(1:p),(p+2)] <- cov_Se2_par
      cov_par[(p+3),c(1:p)] <- cov_par[c(1:p),(p+3)] <- cov_Sp1_par
      cov_par[(p+4),c(1:p)] <- cov_par[c(1:p),(p+4)] <- cov_Sp2_par
      
      cov_par[(p+1),(p+2)] <- cov_par[(p+2),(p+1)] <- cov_Se1_Se2
      cov_par[(p+1),(p+3)] <- cov_par[(p+3),(p+1)] <- cov_Se1_Sp1
      cov_par[(p+1),(p+4)] <- cov_par[(p+4),(p+1)] <- cov_Se1_Sp2
      
      cov_par[(p+2),(p+3)] <- cov_par[(p+3),(p+2)] <- cov_Se2_Sp1
      cov_par[(p+2),(p+4)] <- cov_par[(p+4),(p+2)] <- cov_Se2_Sp2
      
      cov_par[(p+3),(p+4)] <- cov_par[(p+4),(p+3)] <- cov_Sp1_Sp2
      
      cov_VB_par <- solve(-Hessian_manu - cov_par)
      
      VB.mean <- c(results$mu_beta, results$Se, results$Sp)
      #if (sum(is.na(sqrt(diag(cov_VB_par)))) > 0) {
        #print(paste0('Seed ',seed, ": NaNs produced"))}
      VB.std <- c(sqrt(diag(results$var_beta)),
                  (results$Se_alp_beta[1]*results$Se_alp_beta[2])/((results$Se_alp_beta[1]+results$Se_alp_beta[2])^2*(results$Se_alp_beta[1]+results$Se_alp_beta[2]+1)),
                  (results$Se_alp_beta[3]*results$Se_alp_beta[4])/((results$Se_alp_beta[3]+results$Se_alp_beta[4])^2*(results$Se_alp_beta[3]+results$Se_alp_beta[4]+1)),
                  (results$Sp_alp_beta[1]*results$Sp_alp_beta[2])/((results$Sp_alp_beta[1]+results$Sp_alp_beta[2])^2*(results$Sp_alp_beta[1]+results$Sp_alp_beta[2]+1)),
                  (results$Sp_alp_beta[3]*results$Sp_alp_beta[4])/((results$Sp_alp_beta[3]+results$Sp_alp_beta[4])^2*(results$Sp_alp_beta[3]+results$Sp_alp_beta[4]+1)))
      VB.lower <- VB.mean - 1.96*VB.std
      VB.upper <- VB.mean + 1.96*VB.std
      VB.std_cor <- sqrt(abs(diag(cov_VB_par))) 
      VB.lower_cor <- VB.mean - 1.96*VB.std_cor
      VB.upper_cor <- VB.mean + 1.96*VB.std_cor
      VB.alt.low <- c(rep(0,p),
                      qbeta(0.025,results$Se_alp_beta[1],results$Se_alp_beta[2]),
                      qbeta(0.025,results$Se_alp_beta[3],results$Se_alp_beta[4]),
                      qbeta(0.025,results$Sp_alp_beta[1],results$Sp_alp_beta[2]),
                      qbeta(0.025,results$Sp_alp_beta[3],results$Sp_alp_beta[4])
      )
      VB.alt.upp <- c(rep(0,p),
                      qbeta(0.975,results$Se_alp_beta[1],results$Se_alp_beta[2]),
                      qbeta(0.975,results$Se_alp_beta[3],results$Se_alp_beta[4]),
                      qbeta(0.975,results$Sp_alp_beta[1],results$Sp_alp_beta[2]),
                      qbeta(0.975,results$Sp_alp_beta[3],results$Sp_alp_beta[4])
      )
      time2  <- Sys.time()
      VB.output <- data.frame(
                              est_mean = VB.mean,
                              est_std = VB.std,
                              est_lci = VB.lower,
                              est_uci = VB.upper,
                              est_std_cor = VB.std_cor,
                              est_lci_cor = VB.lower_cor,
                              est_uci_cor = VB.upper_cor,
                              VB.alt.low,
                              VB.alt.upp,
                              time = rep(time2 - time1,length(VB.mean)))
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
