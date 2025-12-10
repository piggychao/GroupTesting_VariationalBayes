library("parallel")
library("groupTesting")
library("BayesLogit")
library("MASS")
library("coda")
library("Rcpp")

rm(list=ls())
### Call Rcpp function for sampling y tilde
Rcpp::sourceCpp('Code_MPT.cpp')

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

n_sim <- 500

simu_func <- function(seed){
  set.seed(seed)
  
  ## Part I. Group testing data generation
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
  
  # Part II. Gibbs Sampling
  time1 <- Sys.time()
  gibbs <- function(input_z, input_x, 
                    beta_init, tau_sq, 
                    Se1, Sp1,
                    iter_warmup, iter_sampling){
    
    n <- dim(input_x)[1]
    p <- dim(input_x)[2]
    pool_num <- length(input_z)
    size <- ceiling(n/pool_num)
    
    # pool_ind <- rep(1:pool_num, each=size) #pool index for each individual - only true under current assayID
    pool_ind_cpp <- rep(0:pool_num, each=size) #pool index for each individual, only works for current assayID
    
    prior_var_inv <- (1/tau_sq)*diag(p) #prior var of beta
    prior_mu <- rep(0, p) #prior mu of beta
    
    Se_p <- Se1 #initial pool Se
    Sp_p <- Sp1 #initial pool Sp
    
    #initial z tilde + retest individuals tilde
    # z_til <- rep(0, pool_num)
    # y_til_rt <- rep(0, length(retest_ind))
    
    #initial y tilde
    p_init <- plogis(input_x %*% beta_init)
    y_init <- rbinom(n, 1, p_init)
    
    #results: beta, pool and individual Se&Sp
    results <- matrix(0, nrow = iter_warmup + iter_sampling, ncol = p)
    
    for (g in 1: (iter_warmup + iter_sampling)){
      p_est <- plogis(input_x %*% beta_init)
      
      ## Sampling individual y_i under MPT
      y_tilde <- sampletildey1(pool_ind_cpp, p_est, input_z, y_init, Se_p, Sp_p, size)
      y_init <- y_tilde
      
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
    }
    return(results[(iter_warmup+1):(iter_warmup+iter_sampling),])
  }
  
  # Get estimated parameters
  output <- gibbs(data_z, data_x, 
                  rep(0,dim(data_x)[2]), 1, 
                  Se_true, Sp_true, 
                  2500, 5000)
  time2 <- Sys.time()
  
  mcmcobj <- mcmc(output)
  est_beta_mean <- summary(mcmcobj)$statistics[,"Mean"]
  est_beta_std <- apply(mcmcobj,2, sd)
  est_beta_lci <- summary(mcmcobj)$quantiles[,1]
  est_beta_uci <- summary(mcmcobj)$quantiles[,5]
  
  parameter_true <- param
  output.mcmc <- data.frame(bias = est_beta_mean - parameter_true,
                            std = est_beta_std,
                            coverage = parameter_true >= est_beta_lci & 
                            parameter_true <= est_beta_uci,
                            seed = seed,
                            time = as.numeric(difftime(time2, time1, units = "secs")))
  return(output.mcmc)
}

output <- mclapply(((nm_set-1)*n_sim+1):(nm_set*n_sim), 
                   simu_func, mc.cores = nm_cores)
time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

saveRDS(output, file = paste0('MCMC_MPT_Known_', time,'.rds'))

