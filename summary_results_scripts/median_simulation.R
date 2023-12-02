
# median of the simulations (based on the median number of susceptibes at t=360)

# "spread_summary.txt" for each combination should be generated beforehand 
# using the "summary_results.R" script
######################################################################################

# Load libraries
library(dplyr)
######################################################################################
# Function to find the median simulation based on susceptibles at t=360

sim_median <- function(path){
  df <- read.table(paste0(path,"/spread_summary.txt"))
  dfS_final <- df %>% 
    filter(stat=="Susceptible", time==360) %>%  
    group_by(sim)
  
  sim_summary <- 
    df %>% 
    filter(stat=="Susceptible", time==360) %>%  
    group_by(sim) %>% 
    summarise(final = p) %>% 
    summarise(median(final))
  
  dfS_final$min_dif <- (dfS_final$p-sim_summary$`median(final)`)^2
  sim_m <- dfS_final$sim[dfS_final$min_dif == min(dfS_final$min_dif)]
  spread_m <- read.table(paste0(path,"/spread",sim_m,".txt"))
  survey_m <- read.table(paste0(path,"/surveys",sim_m,".txt"))
  survey_m <- survey_m %>% 
    replace(is.na(.), -5)
  write.txt(spread, paste0(path, "/spread_median.txt"))
  write.txt(survey, paste0(path, "/surveys_median.txt"))
}

###########################################################################
# Apply to each combination of management plans

# Define intervention parameters
surv <- c("one_step_00", "one_step_01", "two_step_00", "two_step_01")
zt <- c(2500, 5000)
er <- c(50, 100)
b <- c(0, 50, 90)

# Generate all combinations of parameters
parameter_combinations <- expand.grid(surv = surv, zt = zt, er = er, b = b)

# Define function to apply for each combination
apply_function <- function(combination) {
  path <- paste0("./results/", combination$surv, "/zt", combination$zt, "_e", combination$er, "_b", combination$b)
  sim_median(path)
}

# Apply the function to each combination
result <- lapply(1:nrow(parameter_combinations), 
                 function(i) apply_function(parameter_combinations[i, ]))

#######################################################################

# Survey effort (number of samples) and inspection intensity (number of cells) 
# of the median of simulations

cells_sampled <- function(path, surv_d, CL, df_0){
  xy_cells <- read.csv("coords.csv")
  surv <- read.table(paste0(path, "/surveys_median.txt"))
  surv$cell <- xy_cells$cell_ha
  
  nc <- ns <- c()
  for(i in 1:(length(surv)-1)){
    nc[i] <- length(unique(surv$cell[surv[,i]>=0]))
    ns[i] <- length(surv[,i][surv[,i]>=0])
  }
  nc <- nc[-1]
  ns <- ns[-1]
  df <- data.frame(time = seq(12, 30*12, 12),
                   n_samples = ns,
                   n_cells = nc)
  
  df$surv <- surv_d
  df$CL <- CL
  df$zt <- str_sub(path, 29, 32)
  df$er <- str_remove(str_sub(path, 35, 37), "_")
  df$b <- str_remove(str_sub(path, -2), "b")
  df_0 <- rbind(df_0,df)
  return(df_0)
}


path <- paste0("./results/", parameter_combinations$surv,
               "/zt", parameter_combinations$zt,
               "_e", parameter_combinations$er,
               "_b", parameter_combinations$b)

surv_d <- rep(surv_d, each=12)
CL <- rep(c("0.8/0.9", "0.99/0.99", "0.8/0.9", "0.6/0.7"), each=12)
df <- data.frame()

# data frame with all
for(i in 1:length(path)){
  df <- cells_sampled(path[i], surv_d[i], CL[i], df)
}

# Write the results to a file
write.table(df, "./results/samples_cells.txt")
