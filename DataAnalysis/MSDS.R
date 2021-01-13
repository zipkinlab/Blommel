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

#Gamma0
mu_s ~ dunif(0, 8)            #Mean
tau_s <- 1/(sig_s * sig_s)    #Precision
sig_s ~ dunif(0, 8)           #Variance

#Gamma1
gamma1 ~ dnorm(0, 0.1)

#Alpha0
mu_a0 ~ dnorm(0, 0.1)        #Mean
tau_a0 ~ dgamma(0.1, 0.1)     #Precision
sig_a0 <- 1/sqrt(tau_a0)      #Variance

#Alpha1
mu_a1 ~ dnorm(0, 0.1)        #Mean
tau_a1 ~ dgamma(0.1, 0.1)     #Precision
sig_a1 <- 1/sqrt(tau_a1)      #Variance

#Alpha2
mu_a2 ~ dnorm(0, 0.1)        #Mean
tau_a2 ~ dgamma(0.1, 0.1)     #Precision
sig_a2 <- 1/sqrt(tau_a2)      #Variance

#Overdispersion
# r.N ~ dunif(0,10)            #Number of groups
omega ~ dunif(0, 1)

for(s in 1:nspec){

#Psi
# tau_p[s] ~ dgamma(0.1, 0.1)  #Precision
# sig_p[s] <- 1/sqrt(tau_p[s]) #Variance

#Sigma
gamma0[s] ~ dnorm(mu_s, tau_s)  #Intercept parameter

#Expected Number of Groups
alpha0[s] ~ dnorm(mu_a0, tau_a0)    #Intercept parameter
alpha1[s] ~ dnorm(mu_a1, tau_a1)    #Effect of region
alpah2[s] ~ dnorm(mu_a2, tau_a2)    #Effect of migration

#------------#
#-LIKELIHOOD-#
#------------#

#-------------------#
#-Distance sampling-#
#-------------------#

for(j in 1:nsites[1]){

# psi[j,s] ~ dnorm(0, tau_p[s])       #Transect effect parameter

#Scale parameter
sigma[j,s] <- exp(gamma0[s] + gamma1 * region[j])

#Construct cell probabilities for nG cells using numerical integration
#Sum of the area (rectangles) under the detection function

for(k in 1:nG){

#Half normal detection function at midpt (length of rectangle)
g[k,j,s] <- exp(-mdpt[k]*mdpt[k]/(2*sigma[j,s]*sigma[j,s]))

#Detection probability for each distance class k (area of each rectangle)
f[k,j,s] <- g[k,j,s] * v/B

#Conditional detection probability (scale to 1)
fc[k,j,s] <- f[k,j,s]/pcap[j,s]

}#end k loop

#Detection probability at each transect (sum of rectangles)
pcap[j,s] <- sum(f[1:nG,j,s])

for(t in nstart[j]:nend[j]){

#Observed population @ each t,j,s (N-mixture)
y[t,j,s] ~ dbin(pcap[j,s], N[t,j,s])

#Latent Number of Groups @ each t,j,s (negative binomial)
N[t,j,s] ~ dpois(lambda.star[t,j,s])

#Expected Number of Groups
# lambda.star[t,j,s] <- rho[t,j,s] * lambda[t,j,s]
lambda.star[t,j,s] <- z[t,j,s] * lambda[t,j,s]
z[t,j,s] ~ dbern(omega)

#Overdispersion parameter for Expected Number of Groups
# rho[t,j,s] ~ dgamma(r.N, r.N)

#Linear predictor for Expected Number of Groups
lambda[t,j,s] <- exp(alpha0[s] + alpha1[s] * region[j] + alpha2[s] * migration[t] + log(offset[j])) # + psi[j,s])

}#end t loop distance sampling

}#end j loop distance sampling

#-----------------#
#-Transect counts-#
#-----------------#

for(j in (nsites[1] + 1):(nsites[1] + nsites[2])){

# psi[j,s] ~ dnorm(0, tau_p[s])       #Transect effect parameter

#Scale parameter
sigma.new[j,s] <- exp(gamma0[s] + gamma1 * region[j])

for(k in 1:8){

#Half normal detection function at midpt (length of rectangle)
g[k,j,s] <- exp(-mdpt[k]*mdpt[k]/(2*sigma.new[j,s]*sigma.new[j,s]))

#Detection probability for each distance class k (area of each rectangle)
f[k,j,s] <- g[k,j,s] * v/B

}#end k loop

#Detection probability at each transect (sum of rectangles)
pdet[j,s] <- sum(f[1:8,j,s])

for(t in nstart[j]:nend[j]){

#Observed population @ each t,j,s (Transect counts)
y[t,j,s] ~ dbin(pdet[j,s], N[t,j,s])

#Latent Number of Groups @ each t,j,s (negative binomial)
N[t,j,s] ~ dpois(lambda.star[t,j,s])

#Expected Number of Groups
# lambda.star[t,j,s] <- rho[t,j,s] * lambda[t,j,s]
lambda.star[t,j,s] <- z[t,j,s] * lambda[t,j,s]
z[t,j,s] ~ dbern(omega)

#Overdispersion parameter for Expected Number of Groups
# rho[t,j,s] ~ dgamma(r.N, r.N)

#Linear predictor for Expected Number of Groups
lambda[t,j,s] <- exp(alpha0[s] + alpha1[s] * region[j] + alpha2[s] * migration[t] + log(offset[j])) # + psi[j,s])

}#end t loop transect counts

}#end j loop transect counts

}#end s loop

for(i in 1:nobs){

#Observed distance classes
dclass[i] ~ dcat(fc[1:nG, site[i], spec[i]])

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

data <- list(y = y, dclass = dclass)

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

inits <- function(){list(mu_s = runif(1, 5, 6), sig_s = runif(1, 0, 1),
                         gamma0 = runif(nspec, 4.75, 6), gamma1 = runif(1, -1, 1),
                         mu_a0 = runif(1, 1, 2), tau_a0 = runif(1, 0, 1), alpha0 = alpha0(),
                         mu_a1 = runif(1, -1, 1), tau_a1 = runif(1, 0, 1), alpha1 = runif(nspec, -1, 1),
                         mu_a2 = runif(1, -1, 1), tau_a2 = runif(1, 0, 1), alpha2 = runif(nspec, -1, 1),
                         # r.N = runif(1, 1, 2), tau_p = runif(nspec, 0, 10),
                         omega = runif(1, 0, 1),
                         N = Nst)}

#--------------------#
#-Parameters to save-#
#--------------------#

params <- c('mu_s', 'sig_s', 'gamma0', 'gamma1',
            'mu_a0', 'sig_a0', 'alpha0', 
            'mu_a1', 'sig_a1', 'alpha1',
            'mu_a2', 'sig_a2', 'alpha2',
            'omega')
            # 'r.N', 'tau_p')

#---------------#
#-MCMC settings-#
#---------------#

model <- nimbleModel(model.code, constants = constants, data = data, inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(c("alpha0", "alpha1", "alpha2"))

MCMCconf$addSampler(target = c("alpha0[1]", "alpha0[2]", "alpha0[3]", "alpha0[4]", "alpha0[5]",
                               "alpha0[6]", "alpha0[7]", "alpha0[8]", "alpha0[9]", "alpha0[10]",
                               "alpha0[11]", "alpha0[12]", "alpha0[13]", "alpha0[14]"),
                    type = "RW_block")

MCMCconf$addSampler(target = c("alpha1[1]", "alpha1[2]", "alpha1[3]", "alpha1[4]", "alpha1[5]",
                               "alpha1[6]", "alpha1[7]", "alpha1[8]", "alpha1[9]", "alpha1[10]",
                               "alpha1[11]", "alpha1[12]", "alpha1[13]", "alpha1[14]"),
                    type = "RW_block")

MCMCconf$addSampler(target = c("alpha2[1]", "alpha2[2]", "alpha2[3]", "alpha2[4]", "alpha2[5]",
                               "alpha2[6]", "alpha2[7]", "alpha2[8]", "alpha2[9]", "alpha2[10]",
                               "alpha2[11]", "alpha2[12]", "alpha2[13]", "alpha2[14]"),
                    type = "RW_block")

nFun <- function(node, s){
  paste(node, "[", s, "]", sep = "")
}

for(s in 1:nspec){
  MCMCconf$addSampler(target = c(nFun("alpha0", s), nFun("alpha1", s), nFun("alpha2", s)),
                      type = "AF_slice")
}

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

ID <- paste("chain", length(list.files(pattern = "chain", full.names = FALSE)) + 1, sep="")
assign(ID, out)
save(list = ID, file = paste0(ID, ".Rds"))
