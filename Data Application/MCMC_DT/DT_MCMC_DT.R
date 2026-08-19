library("parallel")
library("groupTesting")
library("BayesLogit")
library("MASS")
library("coda")
library("Rcpp")
set.seed(2)

rm(list=ls())
### Call Rcpp function for sampling y tilde
Rcpp::sourceCpp('DT_MCMC_DT.cpp')

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

## Dorfman Pool Testing (DT) - Unknown Se&Sp
S <- 2 #2-stage hierarchical testing
psz <- c(psz,1) #Pool sizes used in stages 1-2
Se_true <- c(0.95,0.98) #Sensitivity in stages 1-2
Sp_true <- c(0.98,0.99) #Specificity in stages 1-2
assayID <- c(1,1) #Assays used in stages 1-2

total.mcmc <- 15000
burn.in.step <- 1000
tm.ct <- 1

gibbs <- function(input_z, input_x, input_retest_y, retest_ind,
                  beta_init, tau_sq, 
                  Se1, Se2, Sp1, Sp2, accuracy_init, 
                  iter_warmup, iter_sampling){
    
  n <- dim(input_x)[1]
  p <- dim(input_x)[2]
  pool_num <- length(input_z)
  size <- ceiling(n/pool_num)
  retest_pool_idx <- which(input_z > 0) #positive pool index
  pool_ind_cpp <- rep(0:pool_num, each=size) #pool index for each individual, only works for current assayID
    
  ## Under DT, retest individuals in the positive pools
  y_test <- input_retest_y
  y_test_sub <- y_test[retest_ind]
  
  prior_var_inv <- (1/tau_sq)*diag(p) #prior var of beta
  prior_mu <- rep(0, p) #prior mu of beta
  
  Se_p <- Se1 #initial pool Se
  Se_i <- Se2 #initial individual Se
  Sp_p <- Sp1 #initial pool Sp
  Sp_i <- Sp2 #initial individual Sp
  
  #initial z tilde + retest individuals tilde
  z_til <- rep(0, pool_num)
  y_til_rt <- rep(0, length(retest_ind))
  
  #initial y tilde
  p_init <- plogis(input_x %*% beta_init)
  y_init <- rbinom(n, 1, p_init)
  
  #results: beta, pool and individual Se&Sp
  results <- matrix(0, nrow = iter_warmup + iter_sampling, ncol = (p+4))
  time.vec <- rep(NA,(total.mcmc/burn.in.step + 1))
  time.vec[1] <- Sys.time()
  
  for (g in 1: (iter_warmup + iter_sampling)){
    p_est <- plogis(input_x %*% beta_init)
    
    ## Sampling individual y_i
    y_tilde <- sampletildey2(pool_ind_cpp, p_est, input_z,
                                y_init, input_retest_y,
                                c(Se_p, Se_i), c(Sp_p, Sp_i))
    y_init <- y_tilde
      ## Sampling Sensitivity and Specificity
    ### Pool Se and Sp
    z_matrix <- matrix(y_init, ncol = size, byrow = TRUE) #only works for current assayID
    z_til <- ifelse(rowSums(z_matrix)>0, 1, 0)
      
    Se_p <- rbeta(1, sum(input_z * z_til) + accuracy_init, 
                  sum((1-input_z) * z_til) + accuracy_init) 
    Sp_p <- rbeta(1, sum((1-input_z) * (1-z_til)) + accuracy_init, 
                  sum(input_z * (1-z_til)) + accuracy_init)
    
    ### Individual Se and Sp
    #find the y_tilde for re-test individuals
    y_til_rt <- y_init[retest_ind]
   
    Se_i <- rbeta(1, sum(y_test_sub * y_til_rt) + accuracy_init,
                    sum((1-y_test_sub) * y_til_rt) + accuracy_init)
    Sp_i <- rbeta(1, sum((1-y_test_sub) * (1-y_til_rt)) + accuracy_init,
                    sum(y_test_sub * (1-y_til_rt)) + accuracy_init)
    
    ## Sampling beta
    b_pg <- as.numeric(rep(1,n))
    c_pg <- as.numeric(input_x %*% beta_init)
    omega <- rpg(n, b_pg, c_pg)
    kappa <- y_init - b_pg/2
    var_ome <- prior_var_inv + crossprod(input_x*omega, input_x)
    m_ome <- crossprod(input_x, kappa)
    post_var <- solve(var_ome)
    post_mu <- post_var %*% (m_ome + prior_mu)
    beta_par <- mvrnorm(1, post_mu, Sigma = post_var)
    beta_init <- beta_par
    
    ## save the results
    results[g, 1:p] <- beta_par
    results[g, (p+1):(p+4)] <- c(Se_p, Sp_p, Se_i, Sp_i)
    if(g%%burn.in.step==0){
      print(g)
      tm.ct = tm.ct + 1
      time.vec[tm.ct] = Sys.time()
    }
  }
  return(list(results,time.vec))
}
  
alpha0 <- c(rnorm(1,703,10),rnorm(1,763,10),rnorm(1,3197,10),rnorm(1,1109,10)) #alpha0 for Se_p,Se_i,Sp_p,Sp_i
beta0 <- c(rnorm(1,39,10),rnorm(1,22,10),rnorm(1,65,10),rnorm(1,13,10)) #beta0 for Se_p,Se_i,Sp_p,Sp_i
beta0[beta0<0] <- 0.1
SeSp0 <- alpha0/(alpha0 + beta0) #Initial value for Se & Sp

# Get estimated parameters
output <- gibbs(data_z, data_x, retest_y, retest_ind, 
                rep(0,dim(data_x)[2]), 100, 
                SeSp0[1], SeSp0[2], SeSp0[3], SeSp0[4], 1, 
                0, total.mcmc)

out.mat <- NULL
for(j in 1:(total.mcmc/burn.in.step-1)){
  max.ss = total.mcmc - j*burn.in.step
  f.stop = max.ss/burn.in.step
  for(i in 1:f.stop){
    mcmcobj <- mcmc(output[[1]][(j*burn.in.step+1):(j*burn.in.step+i*burn.in.step),])
    timeobj <- output[[2]]
    est_beta_mean <- summary(mcmcobj)$statistics[,"Mean"]
    est_beta_std <- apply(mcmcobj,2, sd)
    est_beta_lci <- summary(mcmcobj)$quantiles[,1]
    est_beta_uci <- summary(mcmcobj)$quantiles[,5]
    output.mcmc <- rbind(est_beta_mean,
                         est_beta_std,
                         est_beta_lci,
                         est_beta_uci,
                         (j*burn.in.step+1), (i*burn.in.step),
                         timeobj[(j+i)]-timeobj[1], effectiveSize(mcmcobj),  pnorm(abs(geweke.diag(mcmcobj)$z),lower.tail = F))
    out.mat<- rbind(out.mat,output.mcmc)
  }
}


time.mat <- out.mat[seq(from = 7,to = dim(out.mat)[1], by = 9),]


est.mat <- round(out.mat[seq(from = 1,to = dim(out.mat)[1], by = 9),],3)
sd.mat <- round(out.mat[seq(from = 2,to = dim(out.mat)[1], by = 9),],3)
lci.mat <- round(out.mat[seq(from = 3,to = dim(out.mat)[1], by = 9),],3)
uci.mat <- round(out.mat[seq(from = 4,to = dim(out.mat)[1], by = 9),],3)
