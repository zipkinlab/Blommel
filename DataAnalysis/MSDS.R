#------------------------------------------------#
#----Hierarchical Community Distance Sampling----#
#------------------------------------------------#

#----------------#
#-Load Libraries-#
#----------------#

library(nimble)
library(coda)

#-----------#
#-Load Data-#
#-----------#

load(file = "../DataFormatting/FormattedData.Rdata")

#--------------#
#-NIMBLE model-#
#--------------#

model.code <- nimbleCode({
  
  #--------#
  #-PRIORS-#
  #--------#
  
  #Overdispersion
  r.N ~ dunif(0,10) 
  
  #Psi
  # tau_p[s] ~ dgamma(0.1, 0.1)  #Precision
  # sig_p[s] <- 1/sqrt(tau_p[s]) #Variance
  
  #Sigma
  gamma0 ~ dnorm(0, 0.01)  #Intercept parameter
  gamma1 ~ dnorm(0, 0.01)
  
  #Expected Number of Groups
  alpha0 ~ dnorm(0, 0.01)    #Intercept parameter
  alpha1 ~ dnorm(0, 0.01)    #Effect of region
  alpha2 ~ dnorm(0, 0.01)    #Effect of migration
  
  #------------#
  #-LIKELIHOOD-#
  #------------#
  
  #-------------------#
  #-Distance sampling-#
  #-------------------#
  
  for(j in 1:nsites[1]){
    
    # psi[j,s] ~ dnorm(0, tau_p[s])       #Transect effect parameter
    
    #Scale parameter
    sigma[j] <- exp(gamma0 + gamma1 * region[j])
    
    #Construct cell probabilities for nG cells using numerical integration
    #Sum of the area (rectangles) under the detection function
    
    for(k in 1:nG){
      
      #Half normal detection function at midpt (length of rectangle)
      g[k,j] <- exp(-mdpt[k]*mdpt[k]/(2*sigma[j]*sigma[j]))
      
      #Detection probability for each distance class k (area of each rectangle)
      f[k,j] <- g[k,j] * v/B
      
      #Conditional detection probability (scale to 1)
      fc[k,j] <- f[k,j]/pcap[j]
      
    }#end k loop
    
    #Detection probability at each transect (sum of rectangles)
    pcap[j] <- sum(f[1:nG,j])
    
    for(t in nstart[j]:nend[j]){
      
      #Observed population @ each t,j,s (N-mixture)
      y[t,j] ~ dbin(pcap[j], N[t,j])
      
      #Latent Number of Groups @ each t,j,s (negative binomial)
      N[t,j] ~ dpois(lambda.star[t,j])
      
      #Expected Number of Groups
      lambda.star[t,j] <- rho[t,j] * lambda[t,j]
      
      #Overdispersion parameter for Expected Number of Groups
      rho[t,j] ~ dgamma(r.N, r.N)
      
      #Linear predictor for Expected Number of Groups
      lambda[t,j] <- exp(alpha0 + alpha1 * region[j] + alpha2 * migration[t] + log(offset[j])) # + psi[j,s])
      
    }#end t loop distance sampling
    
  }#end j loop distance sampling
  
  for(i in 1:nobs){
    
    #Observed distance classes
    dclass[i] ~ dcat(fc[1:nG, site[i]])
    
  }#end i loop
  
})

#--------------#
#-Compile data-#
#--------------#

attach(Data)

#MTF: for each species, reset nobs, site, y, dclass

constants <- list(nG = nG, v = v, B = B, mdpt = mdpt, nobs = sum(spec==12),
                  nstart = nstart, nend = nend, nsites = nsites,
                  site = site[spec == 12], offset = offset, region = region,
                  migration = migration)

data <- list(y = y[,,12], dclass = dclass[spec == 12])

#----------------#
#-Initial values-#
#----------------#

#MTF: update for each species

Nst <- y[,,12] + 1

#---------------#
#-Inital values-#
#---------------#

inits <- function(){list(gamma0 = runif(1, 2, 6), gamma1 = runif(1, -1, 1),
                         alpha0 = runif(1, 0, 4), alpha1 = runif(1, -1, 1), alpha2 = runif(1, -1, 1),
                         r.N = runif(1, 1, 2), 
                         N = Nst)}

#--------------------#
#-Parameters to save-#
#--------------------#

params <- c('gamma0', 'gamma1',
            'alpha0', 'alpha1', 'alpha2',
            'r.N')

#---------------#
#-MCMC settings-#
#---------------#

model <- nimbleModel(model.code, constants = constants, data = data, inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(c("alpha0", "alpha1", "alpha2"))

MCMCconf$addSampler(target = c("alpha0", "alpha1", "alpha2"), type = "AF_slice")

MCMCconf$removeSampler(c("gamma0", "gamma1"))

MCMCconf$addSampler(target = c("gamma0", "gamma1"), type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

model.comp <- compileNimble(model, MCMC)

nc <- 1
ni <- 50000
nb <- 40000
nt <- 5

out <- runMCMC(model.comp$MCMC, niter = ni, nburnin = nb, nchains = nc, thin = nt, samplesAsCodaMCMC = TRUE)

#-------------#
#-Save output-#
#-------------#

ID <- paste("spec12_chain", length(list.files(pattern = "spec12_chain", full.names = FALSE)) + 1, sep="")
assign(ID, out)
save(list = ID, file = paste0(ID, ".Rds"))
