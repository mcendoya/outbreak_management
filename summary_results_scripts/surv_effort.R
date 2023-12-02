
# Script to extract survey effort (number of samples) and 
# inspection intensity (number of cells)

#################################################################################

# Read the coordinates data
xy_cells <- read.csv("coords.csv")

# Define a function to calculate survey effort and inspection intensity
samples_set <- function(surv_d, zt, er, b, CL){
  sim <- 50
  df_0 <- data.frame()
  
  # Generate the path based on intervention parameters
  path <- paste0("./results/", surv_d,
                 "/zt", zt, "_e", er, "_b", b,
                 "/zt", zt, "_e", er, "_b", b)
  
  # Loop through simulations
  for (s in 1:sim){
    surv <- read.table(paste0(path, "/surveys", s, ".txt"))
    colnames(surv) <- paste0("t", c(0, seq(12, 30*12, 12)))
    surv$cell <- xy_cells$cell_ha
    
    # Replace NaN values with -5
    for(i in 1:length(surv)){
      surv[[i]][is.nan(surv[[i]])] <- -5
    }
    
    # Calculate the number of unique cells and samples at each time point
    nc <- ns <- c()
    for(i in 1:(length(surv)-1)){
      nc[i] <- length(unique(surv$cell[surv[,i]>=0]))
      ns[i] <- length(surv[,i][surv[,i]>=0])
    }
    
    # Summarize data for each simulation
    nc <- sum(nc)
    ns <- sum(ns)
    df <- data.frame(surv = surv_d,
                     CL = CL,
                     zt = zt,
                     er = er,
                     b = b,
                     sim = s,
                     n_samples = ns,
                     n_cells = nc)
    df_0 <- rbind(df_0,df)
  }
  
  # Write the summary data to a file
  write.table(df_0, paste0("./results/",
                           surv_d, "/zt", zt, "_e", er, "_b", b,
                           "/survey_effort.txt"))
}


# Define intervention parameters
zt <- c(2500, 5000)
er <- c(50, 100)
b <- c(0, 50, 90)

# Create intervention lists using pre-defined parameters
interventions <- list(
  list(surv_d = "one_step_00", zt = zt, er = er, b = b, CL = "0.8/0.9"),
  list(surv_d = "one_step_01", zt = zt, er = er, b = b, CL = "0.99/0.99"),
  list(surv_d = "two_step_00", zt = zt, er = er, b = b, CL = "0.8/0.9"),
  list(surv_d = "two_step_01", zt = zt, er = er, b = b, CL = "0.6/0.7")
)

# Loop through interventions and apply the samples_set function
for (intervention in interventions) {
  samples_set(intervention$surv_d, 
              intervention$zt, 
              intervention$er, 
              intervention$b, 
              intervention$CL)
}