sim_output1 <- readRDS("XXX.rds") #input result RDS files
sim_output2 <- readRDS("XXX.rds")
sim_output3 <- readRDS("XXX.rds")
sim_output4 <- readRDS("XXX.rds")
sim_output5 <- readRDS("XXX.rds")
sim_output6 <- readRDS("XXX.rds")
sim_output_list <- list(sim_output1, sim_output2,
                        sim_output3, sim_output4,
                        sim_output5, sim_output6)

settings <- expand.grid(n = c(3000, 6000, 12000), 
                        prev = c(10, 34)) 
for (s in 1:dim(settings)[1]){
  N <- settings[s,1]
  prev <- settings[s,2]
  if (prev==10){
    param <- c(-4,2,1) #Beta_true: prev=10%
  }  
  if (prev==34){
    param <- c(-2,3,1) #Beta_true: prev=34%
  }
  
  sim_output <- sim_output_list[[s]]
  
  mean_est.VB <- matrix(NA, length(sim_output), length(param))

  std_est.VB1 <- matrix(NA, length(sim_output), length(param))
  std_est.VB2 <- matrix(NA, length(sim_output), length(param))
 
  coverage.VB1 <- matrix(NA, length(sim_output), length(param))
  coverage.VB2 <- matrix(NA, length(sim_output), length(param))
 
  time.VB <- rep(NA, length(sim_output))
  
  for (i in 1:length(sim_output)){
    output <- sim_output[[i]]
    
    mean_est.VB[i,] <- output[,1]
    std_est.VB1[i,] <- output[,2]
    coverage.VB1[i,] <- output[,5]
    std_est.VB2[i,] <- output[,6]
    coverage.VB2[i,] <- output[,9]
    time.VB[i] <- as.numeric(output[1,10], units='secs')
    
    if (i %% 50 ==0) {print(output$seed[1])}
  }
  
  output <- data.frame(Bias.VB=apply(mean_est.VB, 2, mean) - param,
                       SE.VB=apply(mean_est.VB, 2, sd),
                       Std.VB=apply(std_est.VB1, 2, mean),
                       Coverage_rate.VB=apply(coverage.VB1, 2, mean),
                       Std.hyd=apply(std_est.VB2, 2, mean),
                       Coverage_rate.hyd=apply(coverage.VB2, 2, mean),
                       Avg_time.VB=mean(time.VB))
  row.names(output) <- c("Beat0", "Beta1", "Beat2")
  
  filename <- paste0("IT_",N, "_prev",prev,".csv") #Change GT protocol name
  
  write.csv(output, file=paste0("Result/",filename), row.names = TRUE)
}
