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
  
  mean_est.MC <- matrix(NA, length(sim_output), length(param))
  std_est.MC <- matrix(NA, length(sim_output), length(param))
  coverage.MC <- matrix(NA, length(sim_output), length(param))
  time.MC <- rep(NA, length(sim_output))
  
  for (i in 1:length(sim_output)){
    output <- sim_output[[i]]
  
    mean_est.MC[i,] <- output$bias
    std_est.MC[i,] <- output$std
    coverage.MC[i,] <- output$coverage
    time.MC[i] <- as.numeric(output$time[1], units='secs')
    if (i %% 50 ==0) {print(output$seed[1])}
  }
  
  output <- data.frame(Bias.MC=apply(mean_est.MC, 2, mean),
                       Std.MC=apply(std_est.MC, 2, mean),
                       Coverage_rate.MC=apply(coverage.MC, 2, mean),
                       Avg_time.MC=mean(time.MC))
  row.names(output) <- c("Beat0", "Beta1", "Beat2")
  
  filename <- paste0("MPT_Knw_",N, "_prev",prev,".csv") #Change GT protocol name
  
  write.csv(output, file=paste0("Result/",filename), row.names = TRUE)
}

