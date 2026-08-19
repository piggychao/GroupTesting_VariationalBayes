library("parallel")
library("groupTesting")
library("BayesLogit")
library("MASS")
library("coda")
library("Rcpp")

### Call Rcpp function for sampling y tilde
Rcpp::sourceCpp('DA_MCMC_MPT.cpp')

##### Printing current parameters ##### 
psz = 3
data_z <- readRDS("Data_z_pz3.rds")
J = length(data_z)
N = J*psz
data_x <- readRDS("Data_x.rds")
data_x <- data_x[1:N,]
p <- dim(data_x)[2]

pool_index_mtx <- readRDS("pool_index_mtx_pz3.rds")


## Master Pool Testing (MPT) - known Se&Sp
psz <- 3 #pool size
Se_true <- 0.9734 #Pool sensitivity
Sp_true <- 0.999 #Pool specificity

total.mcmc <- 15000
burn.in.step <- 1000
tm.ct <- 1
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
    
  pool_ind_cpp <- rep(0:pool_num, each=size) #pool index for each individual, only works for current assayID
  
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
    print(g)
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
                rep(0,dim(data_x)[2]), 100, 
                Se_true, Sp_true, 
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
    output.mcmc <- rbind(est_beta_mean ,
                                est_beta_std,
                                est_beta_lci,
                                est_beta_uci,
                                (j*burn.in.step+1), (i*burn.in.step),
                                timeobj[(j+i)]-timeobj[1], effectiveSize(mcmcobj), pnorm(abs(geweke.diag(mcmcobj)$z),lower.tail = F))
    out.mat<- rbind(out.mat,output.mcmc)
  }
}

time.mat <- out.mat[seq(from = 7,to = dim(out.mat)[1], by = 9),]
est.mat <- round(out.mat[seq(from = 1,to = dim(out.mat)[1], by = 9),],3)
sd.mat <- round(out.mat[seq(from = 2,to = dim(out.mat)[1], by = 9),],3)
lci.mat <- round(out.mat[seq(from = 3,to = dim(out.mat)[1], by = 9),],3)
uci.mat <- round(out.mat[seq(from = 4,to = dim(out.mat)[1], by = 9),],3)
#save.image('C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Tick Surveillance R01/Statistical Methods/Variational Bayes/Data Application/MCMC_MPT/MCMC_MPT_DA.RData')


