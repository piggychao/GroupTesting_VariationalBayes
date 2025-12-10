library(tidyverse)
## Plot box plots for computation time
## 1. Known Se&Sp MPT (VB vs. MCMC)
settings <- expand.grid(n = c(3000, 6000, 12000), 
                        prev = c(10, 34)) 
filename_VB <- c("VB_MPT_known_xxx.rds", # change to output file name
              "VB_MPT_known_xxx.rds",
              "VB_MPT_known_xxx.rds",
              "VB_MPT_known_xxx.rds",
              "VB_MPT_known_xxx.rds",
              "VB_MPT_known_xxx.rds")
filename_MCMC <- c("MCMC_MPT_Known_xxx.rds",
                   "MCMC_MPT_Known_xxx.rds",
                   "MCMC_MPT_Known_xxx.rds",
                   "MCMC_MPT_Known_xxx.rds",
                   "MCMC_MPT_Known_xxx.rds",
                   "MCMC_MPT_Known_xxx.rds")

time_df <- data.frame()
for (i in 1:nrow(settings)){
  VB_df <- readRDS(paste0("./xxx/", filename_VB[i])) # change to output file path
  MCMC_df <- readRDS(paste0("./xxx/", filename_MCMC[i]))
  
  n_sims <- length(VB_df)
  
  time_VB_df <- data.frame(N=rep(settings[i,1],n_sims), 
                           Prevalence=rep(settings[i,2],n_sims),
                           Type="VB",
                           Time=rep(NA,n_sims))
  time_MCMC_df <- data.frame(N=rep(settings[i,1],n_sims), 
                             Prevalence=rep(settings[i,2],n_sims),
                             Type="MCMC",
                             Time=rep(NA,n_sims))
  for (m in 1:n_sims){
    output1 <- VB_df[[m]]
    output2 <- MCMC_df[[m]]
    time_VB_df$Time[m] <- output1$time.VB[1]
    time_MCMC_df$Time[m] <- output2$time[1]
  }
  time_df <- rbind(time_df, time_VB_df, time_MCMC_df)
}

## 2. Known Se&Sp IT (VB vs. MCMC)
settings <- expand.grid(n = c(3000, 6000, 12000), 
                        prev = c(10, 34)) 
filename_VB <- c("VB_IT_known_xxx.rds",
              "VB_IT_known_xxx.rds",
              "VB_IT_known_xxx.rds",
              "VB_IT_known_xxx.rds",
              "VB_IT_known_xxx.rds",
              "VB_IT_known_xxx.rds")
filename_MCMC <- c("MCMC_IT_Known_xxx.rds",
                   "MCMC_IT_Known_xxx.rds",
                   "MCMC_IT_Known_xxx.rds",
                   "MCMC_IT_Known_xxx.rds",
                   "MCMC_IT_Known_xxx.rds",
                   "MCMC_IT_Known_xxx.rds")

time_df2 <- data.frame()
for (i in 1:nrow(settings)){
  VB_df <- readRDS(paste0("./xxx/", filename_VB[i]))
  MCMC_df <- readRDS(paste0("./xxx/", filename_MCMC[i]))
  
  n_sims <- length(VB_df)
  
  time_VB_df <- data.frame(N=rep(settings[i,1],n_sims), 
                           Prevalence=rep(settings[i,2],n_sims),
                           Type="VB",
                           Time=rep(NA,n_sims))
  time_MCMC_df <- data.frame(N=rep(settings[i,1],n_sims), 
                           Prevalence=rep(settings[i,2],n_sims),
                           Type="MCMC",
                           Time=rep(NA,n_sims))
  for (m in 1:n_sims){
    output1 <- VB_df[[m]]
    output2 <- MCMC_df[[m]]
    time_VB_df$Time[m] <- output1$time.VB[1]
    time_MCMC_df$Time[m] <- output2$time[1]
  }
  time_df2 <- rbind(time_df2, time_VB_df, time_MCMC_df)
}

## 3. Unknown Se&Sp DT (VB vs. MCMC)
settings <- expand.grid(n = c(3000, 6000, 12000), 
                        prev = c(10, 34)) 
filename_VB <- c("VB_DT_unknown_xxx.rds",
                 "VB_DT_unknown_xxx.rds",
                 "VB_DT_unknown_xxx.rds",
                 "VB_DT_unknown_xxx.rds",
                 "VB_DT_unknown_xxx.rds",
                 "VB_DT_unknown_xxx.rds")
filename_MCMC <- c("MCMC_DT_unknown_xxx.rds",
                   "MCMC_DT_unknown_xxx.rds",
                   "MCMC_DT_unknown_xxx.rds",
                   "MCMC_DT_unknown_xxx.rds",
                   "MCMC_DT_unknown_xxx.rds",
                   "MCMC_DT_unknown_xxx.rds")

time_df3 <- data.frame()
for (i in 1:nrow(settings)){
  VB_df <- readRDS(paste0("./xxx/", filename_VB[i]))
  MCMC_df <- readRDS(paste0("./xxx/", filename_MCMC[i]))
  n_sims <- length(VB_df)
  
  time_VB_df <- data.frame(N=rep(settings[i,1],n_sims), 
                           Prevalence=rep(settings[i,2],n_sims),
                           Type="VB",
                           Time=rep(NA,n_sims))
  time_MCMC_df <- data.frame(N=rep(settings[i,1],n_sims), 
                             Prevalence=rep(settings[i,2],n_sims),
                             Type="MCMC",
                             Time=rep(NA,n_sims))
  for (m in 1:n_sims){
    output1 <- VB_df[[m]]
    output2 <- MCMC_df[[m]]
    time_VB_df$Time[m] <- output1$time[1]
    time_MCMC_df$Time[m] <- output2$time[1]
  }
  time_df3 <- rbind(time_df3, time_VB_df, time_MCMC_df)
}

time_df_total <- rbind(time_df, time_df2, time_df3)
time_df_total$Setting <- c(rep("MPT Knw.", nrow(time_df)), 
                           rep("IT Knw.", nrow(time_df2)), 
                           rep("DT Unk.", nrow(time_df3)))
time_df_total$N <- as.factor(time_df_total$N)
time_df_total$Prevalence_t[time_df_total$Prevalence == 10] <- "Prevalence = 10%"
time_df_total$Prevalence_t[time_df_total$Prevalence == 34] <- "Prevalence = 34%"
time_df_total$Setting <- factor(time_df_total$Setting, levels = c("MPT Knw.", "IT Knw.", "DT Unk."))

time_plot <- ggplot(time_df_total, aes(x = N, y = Time/60, color = Type)) +
  geom_boxplot() +
  facet_grid(Setting~Prevalence_t) +
  labs(x = "Sample Size", y = "Computation Time (Mintues)") +
  scale_color_manual(values = c("VB" = "blue", "MCMC" = "red"),
                     name = "Method") +
  theme_grey(base_size = 14)
time_plot

ggsave("Time Plot For VB vs MCMC.png", plot = time_plot, width = 8, height = 8, dpi=300)


