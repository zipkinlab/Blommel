#------------------------------------------------#
#----Hierarchical Community Distance Sampling----#
#----Script created by Matthew Farr--------------#
#------------------------------------------------#

#----------#
#-Set Seed-#
#----------#

set.seed(4567)

#---------------#
#-Load Libaries-#
#---------------#

library(jagsUI)

#-----------#
#-Load Data-#
#-----------#

load("herbdata.R")

#---------------#
#-Attach DSdata-#
#---------------#

attach(DSdata)

#-------------------#
#-Compile BUGS data-#
#-------------------#

data <- list(nD = nD, v = v, site = site[spec < 15], rep = rep[spec < 15], spec = spec[spec < 15],
             y = y, B = B, mdpt = mdpt, nherb = nherb, nobs = 22212, dclass = dclass[spec < 15], nsites = nsites, 
             nreps = nreps, gs = gs[spec < 15], Migrant = Migrant, Hyena = Hyena, Lion = Lion, BBJ = BBJ, offset = offset, region = region)

Nst <- y + 1


alpha0 <- function(){
  alpha0 <- c(runif(1,2,3), runif(1,-0.5,0.5), runif(1,1.5,2.5), runif(1,1.5,2.5), runif(1,-1,0), runif(1,-0.5,0.5), runif(1,-0.5,0.5),
              runif(1,2,3), runif(1,2,3), runif(1,3,4), runif(1,3,4), runif(1,-0.5,0.5), runif(1,2,3), runif(1,3,4))
  return(alpha0)
}

beta0 <- function(){
  beta0 <- c(runif(1,3,4), runif(1,1.5,2.5), runif(1,1.5,2.5), runif(1,1,2), runif(1,0.5,1.5), runif(1,1.5,2.5), runif(1,0,1),
              runif(1,2,3), runif(1,2,3), runif(1,1.5,2.5), runif(1,0,1), runif(1,1,2), runif(1,2,3), runif(1,3,4))
  return(beta0)
}


#---------------#
#-Inital values-#
#---------------#

inits <- function(){list(mu_s = runif(1, 5, 6), sig_s = runif(1, 0, 1),
                         gamma0 = runif(nherb, 4.75, 6), gamma1 = runif(1, 0, 0.5),
                         mu_a0 = runif(1, 1, 2), tau_a0 = runif(1, 0, 1), alpha0 = alpha0(), 
                         beta0 = beta0(), r.N = runif(1, 1, 2), r.G = runif(1, 0.9, 1.1), tau_p = runif(nherb, 0, 10),
                         N = Nst)}

#--------------------#
#-Parameters to save-#
#--------------------#

params <- c('mu_s', 'sig_s', 'gamma0', 'gamma1', 
            'mu_a0', 'sig_a0', 'alpha0', 'alpha1', 
            'mu_a2', 'sig_a2', 'alpha2', 
            'mu_a3', 'sig_a3', 'alpha3', 
            'mu_a4', 'sig_a4', 'alpha4', 
            'beta0',
            'r.N', 'r.G', 'tau_p', 'RegGS')

#---------------#
#-MCMC settings-#
#---------------#

nc <- 3
ni <- 50000
nb <- 40000
nt <- 5
na <- 1000

out <- jagsUI(data = data, inits = inits, parameters.to.save = params, model.file = "MS_lulc.txt", 
              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)

save(out, file = "out.Rdata")

