
# Effectiveness and efficiency calculation

# "spread_summary.txt" for each combination should be generated beforehand 
# using the "summary_results.R" script

########################################################

# Load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# Function to find the median simulation based on susceptibles at t=360
sim_median <- function(df){
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
  return(sim_m[1])
}


# Read data for the baseline scenario (no interventions)
df_spread_ni <- read.table("./results/no_intervention/spread_summary.txt")
df_spread_ni_m <- df_spread_ni %>% 
  filter(sim==sim_median(df_spread_ni))

# Extract values for Infected (In) and Susceptible (Sn) at t=360 in the baseline scenario
In <- df_spread_ni_m$n[df_spread_ni_m$stat=="Infected" & df_spread_ni_m$time==360]
Sn <- df_spread_ni_m$n[df_spread_ni_m$stat=="Susceptible" & df_spread_ni_m$time==360]

# Define intervention parameters
surv <- c("one_step_00", "one_step_01", "two_step_00", "two_step_01")
zt <- c(2500, 5000)
er <- c(50, 100)
b <- c(0, 50, 90)

# Generate all combinations of intervention parameters
parameter_combinations <- expand.grid(surv = surv, zt = zt, er = er, b = b)
path <- paste0("./results/", parameter_combinations$surv,
               "/zt", parameter_combinations$zt,
               "_e", parameter_combinations$er,
               "_b", parameter_combinations$b)

# Read spread summary data for each intervention
list_spread <- lapply(path, function(p) read.table(file.path(p, "spread_summary.txt")))
df_spread <- rbindlist(list_spread)

# Find the simulation where the number of S individuals at 360 months is the median
sim_m <- sapply(list_spread, sim_median)

# Filter data for the median simulation for each intervention
for(i in 1:length(list_spread)){
  list_spread[[i]] <-
    list_spread[[i]][list_spread[[i]]$sim==sim_m[i],]
}

df_spread_m <- rbindlist(list_spread)

# Ranking effectiveness and efficiency
df30 <- df_spread_m %>% 
  filter(time==360) %>% 
  select(!c(time, p, sim)) %>% 
  group_by(surv, CL, zt, er, b) %>%  
  summarise(Sc = n[stat=="Susceptible"],
            Ic = n[stat=="Infected"],
            pI= round(((In - Ic)/In)*100,4),
            pS= round(((Sc - Sn)/Sn)*100,4)) %>% 
  ungroup() %>% 
  mutate(rS = dense_rank(-pS),
         rI = dense_rank(-pI)) %>% 
  arrange(rI, rS) %>%
  mutate(Effectiveness = paste0("(", rI, ") ", pI),
         Efficiency = paste0("(", rS, ") ", pS)) %>% 
  select(!c(Sc, Ic, pS, pI, rS, rI)) %>% 
  relocate(surv, CL, zt, er, b, Effectiveness, Efficiency)

# Write the results to a file
write.table(df30, "./results/effectiveness_efficiency.txt", row.names = FALSE)
