#-----------#
#-Libraries-#
#-----------#

library(coda)
library(tidyverse)
library(MCMCvis)

#-----------#
#-Load data-#
#-----------#

#setwd("Z:/Blommel")

pattern <- "spec10_chain"

files <- list.files(path = "./DataAnalysis", pattern = pattern, full.names = TRUE)
#files <- list.files(path = "~/Blommel/DataAnalysis", pattern = pattern, full.names = TRUE)

nc <- length(files)

for(i in 1:nc){
  load(files[i])
}

out <- mcmc.list(mget(paste0(pattern, 1:nc)))

#-------------#
#-Convergence-#
#-------------#

params <- attr(out[[1]], "dimnames")[[2]]

Rhat <- gelman.diag(out[c(1:nc)][,params])
if(all(Rhat[[1]][,1] < 1.1)){
  print("Converged")
}else{
  tmp <- as.numeric(which(Rhat[[1]][,1] > 1.1))
  print("Not converged")
  print(params[tmp])
  traceplot(out[c(1:nc)][,params[tmp]])
}

MCMCsummary(out)


