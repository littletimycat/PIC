library(Rcpp)
library(RcppArmadillo)
source("generate_data_gc.R")
sourceCpp("function_library_gc.cpp")
sourceCpp("function_library.cpp")


SLURM_array_task_id = Sys.getenv('SLURM_ARRAY_TASK_ID')

seed = as.numeric(SLURM_array_task_id)

filename = paste("Updategc/CSV/Updategc_part",seed, ".csv", sep=" ")

#names = as.matrix(c("beta_ipw","gamma_ipw","beta_full","gama_full","ese_beta","ese_gamma","ind_beta", "ind_gamma"), bycol=T)

#write(t(names), file=filename, sep = ",", ncol=length(names),append = F)

source("simulation_gc.R")