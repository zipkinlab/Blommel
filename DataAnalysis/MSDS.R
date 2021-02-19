#------------------------------------------------#
#----Hierarchical Community Distance Sampling----#
#------------------------------------------------#

#---------------------------------#
#-Single Species Integrated Model-#
#---------------------------------#

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

  #Gamma1
  gamma1 ~ dnorm(0, 0.1)


  #Overdispersion
  r.N ~ dunif(0,10)            #Number of groups
  omega ~ dunif(0, 1)
  #tau <- dgamma(0.1,0.1) #ZIF overdispersion


  #Psi
  # tau_p[s] ~ dgamma(0.1, 0.1)  #Precision
  # sig_p[s] <- 1/sqrt(tau_p[s]) #Variance

  #Sigma
  gamma0 ~ dnorm(0, 0.01)  #Intercept parameter

  #Expected Number of Groups
  alpha0 ~ dnorm(0, 0.01)    #Intercept parameter
  alpha1 ~ dnorm(0, 0.01)    #Effect of region
  alpah2 ~ dnorm(0, 0.01)    #Effect of migration

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

      #Observed population @ each t,j (N-mixture)
      y[t,j] ~ dbin(pcap[j], N[t,j])

      #Latent Number of Groups @ each t,j
      #Zero inflated negative binomial (poisson-gamma mixture)
      N[t,j] ~ dpois(lambda.star[t,j])

      #Expected Number of Groups
      # lambda.star[t,j,s] <- rho[t,j,s] * lambda[t,j,s]
      lambda.star[t,j] <- z[t,j] * lambda[t,j] * rho[t,j]
      #zero inflated component
      z[t,j] ~ dbern(omega)

      #Overdispersion parameter for Expected Number of Groups
      rho[t,j] ~ dgamma(r.N, r.N)

      #Linear predictor for Expected Number of Groups
      lambda[t,j] <- exp(alpha0 + alpha1 * region[j] + alpha2 * migration[t] + log(offset[j])) #est[j,s] # + psi[j,s])
      #ZIF over dispersion parameter
      #est[t,j] <- dnorm(0,tau)

    }#end t loop distance sampling

  }#end j loop distance sampling

  #-----------------#
  #-Transect counts-#
  #-----------------#

  for(j in (nsites[1] + 1):(nsites[1] + nsites[2])){

    # psi[j,s] ~ dnorm(0, tau_p[s])       #Transect effect parameter

    #Scale parameter
    sigma.new[j] <- exp(gamma0 + gamma1 * region[j])

    for(k in 1:8){

      #Half normal detection function at midpt (length of rectangle)
      g[k,j] <- exp(-mdpt[k]*mdpt[k]/(2*sigma.new[j]*sigma.new[j]))

      #Detection probability for each distance class k (area of each rectangle)
      f[k,j] <- g[k,j] * v/B

    }#end k loop

    #Detection probability at each transect (sum of rectangles)
    pdet[j] <- sum(f[1:8,j])

    for(t in nstart[j]:nend[j]){

      #Observed population @ each t,j,s (Transect counts)
      y[t,j] ~ dbin(pdet[j], N[t,j])

      #Latent Number of Groups @ each t,j,s (negative binomial)
      N[t,j] ~ dpois(lambda.star[t,j])

      #Expected Number of Groups
      # lambda.star[t,j,s] <- rho[t,j,s] * lambda[t,j,s]
      lambda.star[t,j] <- z[t,j] * lambda[t,j]
      z[t,j] ~ dbern(omega)

      #Overdispersion parameter for Expected Number of Groups
      # rho[t,j,s] ~ dgamma(r.N, r.N)

      #Linear predictor for Expected Number of Groups
      lambda[t,j] <- exp(alpha0 + alpha1 * region[j] + alpha2 * migration[t] + log(offset[j])) # + psi[j,s])

    }#end t loop transect counts

  }#end j loop transect counts



  for(i in 1:nobs){

    #Observed distance classes
    dclass[i] ~ dcat(fc[1:nG, site[i]])

  }#end i loop

})

#--------------#
#-Compile data-#
#--------------#

attach(Data)

constants <- list(nG = nG, v = v, B = B, mdpt = mdpt, nobs = nobs,
                  nstart = nstart, nend = nend, nsites = nsites, nspec = nspec,
                  site = site, spec = spec, offset = offset, region = region,
                  migration = migration)

data <- list(nG = nG, v = v, site = site[spec == 8], rep = rep[spec == 8],
             y = y[,,8], B = B, mdpt = mdpt, nobs = 3003, dclass = dclass[spec == 8], nsites = nsites,
             nreps = nreps, gs = gs[spec == 8], offset = offset, NDVI = NDVI)

#data <- list(y = y, dclass = dclass)

#----------------#
#-Initial values-#
#----------------#

Nst <- y + 1


alpha0 <- function(){
  alpha0 <- c(runif(1,2,3), runif(1,-0.5,0.5), runif(1,1.5,2.5), runif(1,1.5,2.5), runif(1,-1,0), runif(1,-0.5,0.5), runif(1,-0.5,0.5),
              runif(1,2,3), runif(1,2,3), runif(1,3,4), runif(1,3,4), runif(1,-0.5,0.5), runif(1,2,3), runif(1,3,4))
  return(alpha0)
}

#---------------#
#-Inital values-#
#---------------#

inits <- function(){list(
  gamma0 = runif(nspec, 4.75, 6), gamma1 = runif(1, -1, 1),
  alpha0 = alpha0(),
  alpha1 = runif(nspec, -1, 1),
  alpha2 = runif(nspec, -1, 1),
  # r.N = runif(1, 1, 2), tau_p = runif(nspec, 0, 10),
  omega = runif(1, 0, 1),
  N = Nst)}

#--------------------#
#-Parameters to save-#
#--------------------#

params <- c('gamma0', 'gamma1', 'alpha0', 'alpha1', 'alpha2', 'omega')
# 'r.N', 'tau_p')

#---------------#
#-MCMC settings-#
#---------------#

model <- nimbleModel(model.code, constants = constants, data = data, inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

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

ID <- paste("Spec8chain", length(list.files(pattern = "Spec8chain", full.names = FALSE)) + 1, sep="")
assign(ID, out)
save(list = ID, file = paste0(ID, ".Rds"))
