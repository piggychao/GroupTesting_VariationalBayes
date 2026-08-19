library("parallel")
library("groupTesting")
library("BayesLogit")
library("MASS")
library("coda")
library("Rcpp")

### Call Rcpp function for sampling y tilde
#Rcpp::sourceCpp('Sims/IT/MCMC/Code_IT.cpp')
Rcpp::sourceCpp('/work/scwatson/VB/IT_MCMC/Code_IT.cpp')

##### Printing current parameters ##### 
N <- 12000 #Sample size

param <- c(-4,2,1) #Beta_true: prev=10%

## Individual Pool Testing (IT) - known Se&Sp
S <- 1 #1-stage hierarchical testing
psz <- 1 #Pool sizes used in stages 1
Se_true <- 0.98 #Pool sensitivity
Sp_true <- 0.99 #Pool specificity
assayID <- 1 #Assays used in individual pools
tau_sq <- 100
n_sim <- 5

total.mcmc <- 15000
burn.in.step <- 1000
tm.ct <- 1

simu_func <- function(seed){
  #set.seed(seed)
  
  ## Part I. Group testing data generation
  x <- cbind(1, rnorm(N), rbinom(N,1,0.5))
  colnames(x) <- c("Intercept", "x1", "x2")
  
  pReg <- plogis(x %*% param)
  # Simulating test responses
  IT.gtData <- hier.gt.simulation(N, p=pReg, S, psz, Se=Se_true, Sp=Sp_true, assayID)$gtData
  
  J <- N/psz[1]
  # pool_index_mtx <- DT.gtData[c(1:J),c(6:8)]
  data_x <- x
  data_z <- IT.gtData[c(1:J),1]
  
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
    
    prior_var_inv <- (1/tau_sq)*diag(p) #prior var of beta
    prior_mu <- rep(0, p) #prior mu of beta
    
    Se_p <- Se1 #initial pool Se
    Sp_p <- Sp1 #initial pool Sp
    
    #initial y tilde
    p_init <- plogis(input_x %*% beta_init)
    y_init <- rbinom(n, 1, p_init)
    
    #results: beta, pool and individual Se&Sp
    results <- matrix(0, nrow = iter_warmup + iter_sampling, ncol = p)
    time.vec <- rep(NA,(total.mcmc/burn.in.step + 1))
    time.vec[1] <- Sys.time()
    for (g in 1: (iter_warmup + iter_sampling)){
      p_est <- plogis(input_x %*% beta_init)
      
      ## Sampling individual y_i
      y_tilde <- sampletildey_it(p_est, input_z, Se_p, Sp_p)
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
      if(g%%burn.in.step==0){
        print(g)
        tm.ct = tm.ct + 1
        time.vec[tm.ct] = Sys.time()
      }
    }
    return(list(results,time.vec))
 }

  # Get estimated parameters
  output <- gibbs(data_z, data_x,  
                  rep(0,dim(data_x)[2]), tau_sq, 
                  Se_true, Sp_true, 
                  0, total.mcmc)
  time2 <- Sys.time()
  
  mcmcobj <- mcmc(output[[1]])
  est_beta_mean <- summary(mcmcobj)$statistics[,"Mean"]
  est_beta_std <- apply(mcmcobj,2, sd)
  est_beta_lci <- summary(mcmcobj)$quantiles[,1]
  est_beta_uci <- summary(mcmcobj)$quantiles[,5]
  convg <- pnorm(abs(geweke.diag(mcmcobj)$z), lower.tail = F)*2 
  
  time2 <- Sys.time()
  
  return(output)
}

out.mat <- NULL
parameter_true <- param
for(sims in 1:n_sim){
  output = simu_func(1)
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
      output.mcmc <- rbind(est_beta_mean - parameter_true,
                           est_beta_std,
                           parameter_true >= est_beta_lci & parameter_true <= est_beta_uci,
                           (j*burn.in.step+1), (i*burn.in.step),
                           timeobj[(j+i)]-timeobj[1], effectiveSize(mcmcobj))
      out.mat<- rbind(out.mat,output.mcmc)
    }
  }
  print(sims)
}

write.table(as.matrix(out.mat),'betares.txt',row.names = FALSE, col.names = FALSE)
