#####################################################################################
# Script to create summary of simulation results
#####################################################################################

# Load libraries
library(stringr)
#####################################################################################

#Baseline scenario (No interventions)

# Set the path for the results folder without interventions
path <- "./results/no_intervention"
df_sims <- data.frame()

# Loop through 50 simulations
for (s in 1:50){
  spread <- read.table(paste0(path, "/spread",s,".txt"))
  colnames(spread) <- paste0("t", seq(0, 30*12, 6))
  
  # Calculate infected and susceptible counts
  indI <- apply(spread, 2, function(x) length(x[x>0])) # Infected (asymptomatic and symptomatic)
  indS <- apply(spread, 2, function(x) length(x[x==0]))# Susceptible
  
  # Create a data frame for each simulation and concatenate to df_sims
  df <- data.frame(time = c(rep(seq(0,30*12,6), 2)),
                   n = c(indI, indS),
                   stat = c(rep("Infected", length(indI)),
                            rep("Susceptible", length(indS))))
  
  df$p <- (df$n/nrow(spread))*100 # percentage
  df$sim <- s
  df_sims <- rbind(df_sims, df)
}

# Write the summary data to a file
write.table(df_sims, paste0(path, "/spread_summary.txt"))

#####################################################################################

# Function to summarize simulation data with interventions

summary_data_simulations <- function(path, surv_d, CL, zt, er, b, sim=50){
  df_sims <- data.frame()
  df_surv_sims <- data.frame()
  df_E_sims <- data.frame()
  
  # Loop through simulations
  for (s in 1:sim){
    spread <- read.table(paste0(path, "/spread", s, ".txt"))
    surv <- read.table(paste0(path, "/surveys", s, ".txt"))
    
    # Process spread and survey data
    colnames(spread) <- paste0("t", seq(0, 30*12, 6))
    colnames(surv) <- paste0("t", c(0, seq(12, 30*12, 12)))
    
    # Replace NaN values in survey data with -5
    for(i in 1:length(surv)){
      surv[[i]][is.nan(surv[[i]])] <- -5
    }
    
    # Calculate infected, susceptible, and eradicated counts
    indI <- apply(spread, 2, function(x) length(x[x>0])) # Infected (asymptomatic and symptomatic)
    indS <- apply(spread, 2, function(x) length(x[x==0]))# Susceptible
    indE <- apply(spread, 2, function(x) length(x[x==-1])) # Eradicated
    
    # Create a data frame for each simulation and concatenate to df_sims
    df <- data.frame(time = c(rep(seq(0,30*12,6), 3)),
                     n = c(indI, indS, indE),
                     stat = c(rep("Infected",length(indI)),
                              rep("Susceptible", length(indS)),
                              rep("Eradicated", length(indE))))
    df$p <- (df$n/nrow(spread))*100 # percentage
    
    # Add intervention parameters and simulation information
    df$surv <- surv_d
    df$CL <- CL
    df$zt <- zt
    df$er <- er
    df$b <- b
    df$sim <- s
    
    df_sims <- rbind(df_sims, df)
    
    # Process survey data
    t_surv <- c(0, seq(12, 30*12, 12))
    indI_s <- apply(surv, 2, function(x) length(x[x>0])) # Infected (asymptomatic and symptomatic) sampled
    indI_ts <- df$n[df$stat=="Infected"&df$time%in%t_surv] # Infected (asymptomatic and symptomatic)
    indS_s <- apply(surv, 2, function(x) length(x[x==0])) # Susceptible sampled
    indS_ts <- df$n[df$stat=="Susceptible"&df$time%in%t_surv]# Susceptible 
    
    # Create a data frame for each simulation and concatenate to df_surv_sims
    df_surv <- data.frame(time = rep(t_surv,2),
                          stat = c(rep("Infected", length(t_surv)),
                                   rep("Susceptible", length(t_surv))),
                          n = c(indI_ts, indS_ts),
                          samples = c(indI_s, indS_s),
                          p_detected = c((indI_s/indI_ts)*100,
                                         (indS_s/indS_ts)*100))
    df_surv$surv <- surv_d
    df_surv$CL <- CL
    df_surv$zt <- zt
    df_surv$er <- er
    df_surv$b <- b
    df_surv$sim <- s
    df_surv_sims <- rbind(df_surv_sims, df_surv)
    
    # Process eradication data
    t <- t_surv[1:length(seq(12, 12*30, by=12))]
    te <- t+6
    se <- c()
    ie <- c()
    
    # Calculate eradicated susceptible and infected counts
    for(i in 1:length(t)){
      t_i <- paste0("t", t[i])
      te_i <- paste0("t", te[i])
      se[i] <- length(spread[,t_i][spread[,te_i] == -1 & spread[,t_i] == 0]) # Eradicated susceptible
      ie[i] <- length(spread[,t_i][spread[,te_i] == -1 & spread[,t_i] > 0]) # Eradicated infected
    }
    
    # Create a data frame for each simulation and concatenate to df_E_sims
    df_E <- data.frame(time = rep(te, 2),
                       n = c(se, ie),
                       stat = c(rep("Susceptible", length(se)),
                                rep("Infected", length(ie))))
    df_E$surv <- surv_d
    df_E$CL <- CL
    df_E$zt <- zt
    df_E$er <- er
    df_E$b <- b
    df_E$sim <- s
    df_E_sims <- rbind(df_E_sims, df_E)
  }
  
  # Write summary data to files
  write.table(df_sims, paste0(path,"/spread_summary.txt"))
  write.table(df_surv_sims, paste0(path,"/survey_summary.txt"))
  write.table(df_E_sims, paste0(path,"/eradication.txt"))
}
#############################################################################################

# Summarize the data for each combination of outbreak management plans.
#############################################################################################

summarize_simulations <- function(strategy, surv_d, CL, zt_values, er_values, b_values) {
  # Loop through the intervention parameters and summarize simulations
  for (r_zt in zt_values) {
    for (r_er in er_values) {
      for (b_i in b_values) {
        # Generate the path
        path <- paste0("./results/", strategy, "/zt", r_zt, "_e", r_er, "_b", b_i)
        # Call the summary_data_simulations function
        summary_data_simulations(path, 
                                 surv_d, 
                                 CL, 
                                 as.character(r_zt), 
                                 as.character(r_er), 
                                 as.character(b_i), 
                                 50) # 50 = number of simulations
      }
    }
  }
}

# Define intervention parameters
zt_values <- c(2500, 5000)
er_values <- c(50, 100)
b_values <- c(0, 50, 90)

# Summarize simulations for different strategies, surv_d, and CL
summarize_simulations("one_step_00", "One-step", "0.8/0.9", zt_values, er_values, b_values)
summarize_simulations("one_step_01", "One-step", "0.99/0.99", zt_values, er_values, b_values)
summarize_simulations("two_step_00", "Two-step", "0.8/0.9", zt_values, er_values, b_values)
summarize_simulations("two_step_01", "Two-step", "0.6/0.7", zt_values, er_values, b_values)





