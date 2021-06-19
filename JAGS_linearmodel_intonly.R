#Kayla Blincow
#Update 3/23/2021

#Fit a linear model in JAGS for C and N v. length

#3/23/2021 Include a random effect of catchsite and year
#4/5/2021 convert single categorical variable models to intercept only models


#clear my workspace
rm(list = ls())

#load necessary libraries
library(R2jags)
library(superdiag)
library(tidyverse)
library(lme4)
library(mcmcplots)
library(viridis)
library(MCMCvis)
library(loo)

#read in our data
d <- read.csv("FinalGSBBulk.csv", header = T)

#convert sample_loc to a factor
d$sample_loc <- factor(d$sample_loc)

#scale data so there is a mean of 0 and a sd of 1
d$Cstd <- as.numeric(scale(d$d13C))
d$Nstd <- as.numeric(scale(d$d15N))


#set MCMC settings for the whole deal
ni <- 350000
nt <- 25
nb <- 25000
nc <- 3

#generate sample index for plotting
smple <- sample(1:39000, 1000)

#create categorical variables
site <- factor(d$sample_loc)

year <- factor(d$year)

d$yr_site <- paste(d$year, d$sample_loc)
yr_site <- factor(d$yr_site)



# Now let's do the models in JAGS! 
#Create the JAGS files

####Iso v. Length####

# CARBON
sink("C_v_Length.txt")
cat("
    model {
    
    # Priors
    
    alpha ~ dnorm(0, 1)
    beta ~ dnorm(0, 1)
    
    sigma ~ dexp(1)
    tau <- 1/(sigma*sigma)
    
    # Likelihood
    
    for(i in 1:n){
    C[i] ~ dnorm(exp.C[i], tau)
    exp.C[i] <- alpha + beta * TL[i]
    LogLik[i] <- logdensity.norm(C[i], exp.C[i], tau)
    }
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(C = d$Cstd, n = length(d$d13C), TL = d$TotalLength)

# Initial values
inits <- function() list(alpha = runif(1, -0.5, .5), 
                         sigma = runif(1, 0, 1), 
                         beta = runif(1, -2, 2))

# Parameters you want to keep track of:
jags.params <- c("alpha", "beta", "sigma", "LogLik")


# Call JAGS from R
mod.C.TL = jags(jags.data,
                inits = NULL,
                parameters.to.save= jags.params,
                model.file="C_v_Length.txt",
                n.chains = nc,
                n.burnin =nb,
                n.thin = nt,
                n.iter = ni,
                DIC = TRUE)

mod.C.TL

#look at parameter estimates
MCMCplot(mod.C.TL,
         params = "alpha")
MCMCplot(mod.C.TL,
         params = "beta")


#PLOT IT
ggplot(data = d, aes(x = TotalLength, y = Cstd)) +
  geom_point() +
  geom_abline(intercept = mod.C.TL$BUGSoutput$sims.list$alpha[smple], 
              slope = mod.C.TL$BUGSoutput$sims.list$beta[smple], 
              color = "gray", alpha = 0.2)+
  geom_abline(slope = mod.C.TL$BUGSoutput$mean$beta, 
              intercept = mod.C.TL$BUGSoutput$mean$alpha, 
              size  = 1, color = "turquoise4") +
  labs(x = "Total Length (cm)", y = expression(delta~" "^13~C)) +
  theme_classic()


#NITROGEN#

sink("N_v_Length.txt")
cat("
    model {
    
    # Priors
    
    alpha ~ dnorm(0, 1)
    beta ~ dnorm(0, 1)
    
    sigma ~ dexp(1)
    tau <- 1/(sigma*sigma)
    
    # Likelihood
    
    for(i in 1:n){
    exp.N[i] <- alpha + beta * TL[i]
    N[i] ~ dnorm(exp.N[i], tau) 
    LogLik[i] <- logdensity.norm(N[i], exp.N[i], tau)
    }
    
    }
    ", fill = TRUE)
sink()



# Bundle data
jags.data <- list(N = d$Nstd, n = length(d$d15N), TL = d$TotalLength)

# Initial values
inits <- function() list(alpha = runif(1, -0.5, .5), 
                         sigma = runif(1, 0, 1), 
                         beta = runif(1, -2, 2))

# Parameters you want to keep track of:
jags.params <- c("alpha", "beta", "sigma", "LogLik")



# Call JAGS from R
mod.N.TL = jags(jags.data,
                inits = NULL,
                parameters.to.save= jags.params,
                model.file="N_v_Length.txt",
                n.chains = nc,
                n.burnin =nb,
                n.thin = nt,
                n.iter = ni,
                DIC = TRUE)

mod.N.TL

#look at parameter estimates
MCMCplot(mod.N.TL,
         params = "alpha")
MCMCplot(mod.N.TL,
         params = "beta")


#PLOT IT
ggplot(data = d, aes(x = TotalLength, y = Nstd)) +
  geom_point() +
  geom_abline(intercept = mod.N.TL$BUGSoutput$sims.list$alpha[smple], 
              slope = mod.N.TL$BUGSoutput$sims.list$beta[smple], 
              color = "gray", alpha = 0.2)+
  geom_abline(slope = mod.N.TL$BUGSoutput$mean$beta, 
              intercept = mod.N.TL$BUGSoutput$mean$alpha, 
              size  = 1, color = "turquoise4") +
  labs(x = "Total Length (cm)", y = expression(delta~" "^13~C)) +
  theme_classic()


####ISO V. SITE####


#CARBON
#set up the JAGS model
sink("C_v_site.txt")
cat("
    model {
    

    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of catch location (i.e. site)
    for(j in 1:Nsite){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    # Likelihood
    for(i in 1:n){
    exp.C[i] <- mu.alpha + alpha[site[i]]
    C[i] ~ dnorm(exp.C[i], tau)  
    LogLik[i] <- logdensity.norm(C[i], exp.C[i], tau)
    }
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(C = d$Cstd, n = length(d$d13C), site = as.numeric(site),
                  Nsite = length(unique(site)))

# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")


# Call JAGS from R
mod.C.site = jags(jags.data,
                  inits = NULL,
                  parameters.to.save= jags.params,
                  model.file="C_v_site.txt",
                  n.chains = nc,
                  n.burnin =nb,
                  n.thin = nt,
                  n.iter = ni,
                  DIC = TRUE)

mod.C.site

#look at parameter estimates
MCMCplot(mod.C.site,
         params = "alpha")


#calculate total variance
tot_var <- mod.C.site$BUGSoutput$mean$sigma^2 + 
  mod.C.site$BUGSoutput$mean$sigma.alpha^2

#calculate variance associated with catch site parameter
site_var_ratio <- (mod.C.site$BUGSoutput$mean$sigma.alpha^2)/tot_var

#NITROGEN#
#set up the JAGS model
sink("N_v_site.txt")
cat("
    model {
    
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of catch location (i.e. site)
    for(j in 1:Nsite){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    # Likelihood
    for(i in 1:n){
    exp.N[i] <- mu.alpha + alpha[site[i]]
    N[i] ~ dnorm(exp.N[i], tau)  
    LogLik[i] <- logdensity.norm(N[i], exp.N[i], tau)
    }
  }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(N = d$Nstd, n = length(d$d15N), site = as.numeric(site),
                  Nsite = length(unique(site)))

# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")

# Call JAGS from R
mod.N.site = jags(jags.data,
                  inits = NULL,
                  parameters.to.save= jags.params,
                  model.file="N_v_site.txt",
                  n.chains = nc,
                  n.burnin =nb,
                  n.thin = nt,
                  n.iter = ni,
                  DIC = TRUE)

mod.N.site

#look at parameter estimates
MCMCplot(mod.N.site,
         params = "alpha")


####ISO V. Year####


#CARBON
#set up the JAGS model
sink("C_v_year.txt")
cat("
    model {
    
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of year
    for(j in 1:Nyear){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    

    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    # Likelihood
    for(i in 1:n){
    exp.C[i] <- mu.alpha + alpha[year[i]]
    C[i] ~ dnorm(exp.C[i], tau)  
    LogLik[i] <- logdensity.norm(C[i], exp.C[i], tau)
    }
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(C = d$Cstd, n = length(d$d13C), year = as.numeric(year),
                  Nyear = length(unique(year)))

# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")



# Call JAGS from R
mod.C.year = jags(jags.data,
                  inits = NULL,
                  parameters.to.save= jags.params,
                  model.file="C_v_year.txt",
                  n.chains = nc,
                  n.burnin =nb,
                  n.thin = nt,
                  n.iter = ni,
                  DIC = TRUE)

mod.C.year

#look at parameter estimates
MCMCplot(mod.C.year,
         params = "alpha")



#NITROGEN
#set up the JAGS model
sink("N_v_year.txt")
cat("
        model {
    
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of year
    for(j in 1:Nyear){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    # Likelihood
    for(i in 1:n){
    exp.N[i] <- mu.alpha + alpha[year[i]]
    N[i] ~ dnorm(exp.N[i], tau)  
    LogLik[i] <- logdensity.norm(N[i], exp.N[i], tau)
    }
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(N = d$Nstd, n = length(d$d15N), year = as.numeric(year),
                  Nyear = length(unique(year)))

# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")


# Call JAGS from R
mod.N.year = jags(jags.data,
                  inits = NULL,
                  parameters.to.save= jags.params,
                  model.file="N_v_year.txt",
                  n.chains = nc,
                  n.burnin =nb,
                  n.thin = nt,
                  n.iter = ni,
                  DIC = TRUE)

mod.N.year

#look at parameter estimates
MCMCplot(mod.N.year,
         params = "alpha")



####ISO V. YearSite####
#CARBON

#set up the JAGS model
sink("C_v_yrsite.txt")
cat("
    model {
    
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of catch location/year
    for(j in 1:Nyrsite){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    # Likelihood
    for(i in 1:n){
    exp.C[i] <- mu.alpha + alpha[yr_site[i]]
    C[i] ~ dnorm(exp.C[i], tau)  
    LogLik[i] <- logdensity.norm(C[i], exp.C[i], tau)
    }
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(C = d$Cstd, n = length(d$d13C), yr_site = as.numeric(yr_site),
                  Nyrsite = length(unique(yr_site)))

# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")

# Call JAGS from R
mod.C.yrsite = jags(jags.data,
                    inits = NULL,
                    parameters.to.save= jags.params,
                    model.file="C_v_yrsite.txt",
                    n.chains = nc,
                    n.burnin =nb,
                    n.thin = nt,
                    n.iter = ni,
                    DIC = TRUE)

mod.C.yrsite

#look at parameter estimates
MCMCplot(mod.C.yrsite,
         params = "alpha")


#NITROGEN
#set up the JAGS model
sink("N_v_yrsite.txt")
cat("
    model {
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of catch location/year
    for(j in 1:Nyrsite){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    # Likelihood
    for(i in 1:n){
    exp.N[i] <- mu.alpha + alpha[yr_site[i]]
    N[i] ~ dnorm(exp.N[i], tau)  
    LogLik[i] <- logdensity.norm(N[i], exp.N[i], tau)

    }
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(N = d$Nstd, n = length(d$d15N), yr_site = as.numeric(yr_site),
                  Nyrsite = length(unique(yr_site)))

# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")

# Call JAGS from R
mod.N.yrsite = jags(jags.data,
                    inits = NULL,
                    parameters.to.save= jags.params,
                    model.file="N_v_yrsite.txt",
                    n.chains = nc,
                    n.burnin =nb,
                    n.thin = nt,
                    n.iter = ni,
                    DIC = TRUE)

mod.N.yrsite

#look at parameter estimates
MCMCplot(mod.N.yrsite,
         params = "alpha")




#####ISO v. Length + (random effect of catch location)####
#CARBON
# Specify model in BUGS language
sink("C_v_Length_site.txt")
cat("
    model {
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of catch location (i.e. site)
    for(j in 1:Nsite){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    beta ~ dnorm(0, 1)

    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    
    # Likelihood
    
    for(i in 1:n){
    C[i] ~ dnorm(exp.C[i], tau)
    exp.C[i] <- mu.alpha + alpha[site[i]] + beta*TL[i] 
    LogLik[i] <- logdensity.norm(C[i], exp.C[i], tau)
    }
    
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(C = d$Cstd, 
                  TL = d$TotalLength,
                  n = length(d$d13C), 
                  site = as.numeric(site),
                  Nsite = length(unique(site)))


# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         beta = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("beta", "alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")


# Call JAGS from R
mod.C.TLsite = jags(jags.data,
                    inits = NULL,
                    parameters.to.save= jags.params,
                    model.file = "C_v_Length_site.txt",
                    n.chains = nc,
                    n.burnin =nb,
                    n.thin = nt,
                    n.iter = ni,
                    DIC = TRUE)

#check model output 
#check for convergence, compare sigma estimates
mod.C.TLsite
#mcmcplot(mod.C.TLsite) #takes some time

#calculate total variance
tot_var <- mod.C.TLsite$BUGSoutput$mean$sigma + 
  mod.C.TLsite$BUGSoutput$mean$sigma.alpha

#calculate variance associated with catch site parameter
site_var_ratio <- (mod.C.TLsite$BUGSoutput$mean$sigma.alpha)/tot_var


#plot my results..
MCMCplot(mod.C.TLsite,
         params = c("alpha", "mu.alpha"))
MCMCplot(mod.C.TLsite,
         params = c("beta"))

#convert JAGS output into new format which is easier to plot

#create objects out of output
mn_site_alpha <- mod.C.TLsite$BUGSoutput$mean$alpha
drws_site_alpha <- mod.C.TLsite$BUGSoutput$sims.list$alpha

mn_mu_alpha <- mod.C.TLsite$BUGSoutput$mean$mu.alpha
drws_mu_alpha <- mod.C.TLsite$BUGSoutput$sims.list$mu.alpha

mn_beta <- mod.C.TLsite$BUGSoutput$mean$beta
drws_beta <- mod.C.TLsite$BUGSoutput$sims.list$beta

#combine alpha estimates into a single mean estimate for each site
alpha_musite <- matrix(NA, nrow = nrow(mn_site_alpha))
for(i in 1:nrow(mn_site_alpha)){
  alpha_musite[i] <- mn_site_alpha[i] + mn_mu_alpha  
}

#combine for the posterior draws as well
alpha_drwssite <- matrix(NA, nrow(drws_site_alpha), 
                         ncol = ncol(drws_site_alpha))
for(i in 1:ncol(drws_site_alpha)){
  alpha_drwssite[,i] <- drws_site_alpha[,i] + drws_mu_alpha  
}

#assign sample sites as column names
colnames(alpha_drwssite) <- c(1:13)

#create a new dataframe for plotting
post_data <- alpha_drwssite %>% 
  as_tibble(.name_repair = "unique") %>% 
  pivot_longer(cols = 1:13, names_to = "site", values_to = "alpha_drws")
post_data$alphamn <- rep(alpha_musite, 39000)
post_data$betadrws <- rep(drws_beta, 13)
post_data$betamn <- rep(mn_beta, nrow(post_data))

#plot my results some more..
colors <- viridis(13)
lsz <- 0.75
smpl <- sample(nrow(post_data), 1000)
d$sample_loc <- factor(d$sample_loc, 
                       labels = c("Ajusco Erendira", "Bahia de Los Angeles",
                                  "Bahia Tortugas", "Carlsbad", "El Datil",
                                  "El Rosario", "Guerrero Negro", 
                                  "Isla Natividad", "Laguna San Ignacio",
                                  "La Jolla", "Punta Abreojos", "San Juanico",
                                  "Solana Beach"))

CvLvS <- ggplot() +
  scale_color_manual(values = colors) +
  geom_abline(intercept = post_data$alpha_drws[smpl],
              slope = post_data$betadrws[smpl],
              color = "gray70", alpha = 0.2) +
  geom_point(data = d, aes(x = TotalLength, y = Cstd, color = sample_loc)) +
  geom_abline(intercept = post_data$alphamn[1],
              slope = post_data$betamn[1], 
              color = colors[1], size = lsz)+
  geom_abline(intercept = post_data$alphamn[2],
              slope = post_data$betamn[2], 
              color = colors[2], size = lsz)+
  geom_abline(intercept = post_data$alphamn[3],
              slope = post_data$betamn[3], 
              color = colors[3], size = lsz)+
  geom_abline(intercept = post_data$alphamn[4],
              slope = post_data$betamn[4], 
              color = colors[4], size = lsz)+
  geom_abline(intercept = post_data$alphamn[5],
              slope = post_data$betamn[5], 
              color = colors[5], size = lsz)+
  geom_abline(intercept = post_data$alphamn[6],
              slope = post_data$betamn[6], 
              color = colors[6], size = lsz)+
  geom_abline(intercept = post_data$alphamn[7],
              slope = post_data$betamn[7], 
              color = colors[7], size = lsz)+
  geom_abline(intercept = post_data$alphamn[8],
              slope = post_data$betamn[8], 
              color = colors[8], size = lsz)+
  geom_abline(intercept = post_data$alphamn[9],
              slope = post_data$betamn[9], 
              color = colors[9], size = lsz)+
  geom_abline(intercept = post_data$alphamn[10],
              slope = post_data$betamn[10], 
              color = colors[10], size = lsz)+
  geom_abline(intercept = post_data$alphamn[11],
              slope = post_data$betamn[11], 
              color = colors[11], size = lsz)+
  geom_abline(intercept = post_data$alphamn[12],
              slope = post_data$betamn[12], 
              color = colors[12], size = lsz)+
  geom_abline(intercept = post_data$alphamn[13],
              slope = post_data$betamn[13], 
              color = colors[13], size = lsz)+
  labs(x = "Total Length (cm)", y = expression(delta~""^13~C),
       color = "Sample Site") +
  theme_classic()

png(filename="MSFigures/CvTLvSite.png", 
    units="in", 
    width=7, 
    height=6, 
    pointsize=8, 
    res=400)
CvLvS

dev.off()

#Or make everyone's lives easier and just plot mu alpha estimate??
pCvL <- ggplot() +
  geom_abline(intercept = drws_mu_alpha[smple],
              slope = drws_beta[smple],
              color = "gray70", alpha = 0.2) +
  geom_point(data = d, aes(x = TotalLength, y = Cstd)) +
  geom_abline(intercept = mn_mu_alpha,
              slope = mn_beta, 
              size = lsz) +
  labs(x = "Total Length (cm)", y = expression(delta~""^13~C)) +
  theme_classic()

png(filename="MSFigures/CvTL.png", 
    units="in", 
    width=6, 
    height=6, 
    pointsize=8, 
    res=400)
pCvL

dev.off()

#NITROGEN
sink("N_v_Length_site.txt")
cat("
    model {
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of catch location (i.e. site)
    for(j in 1:Nsite){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    beta ~ dnorm(0, 1)
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    
    # Likelihood
    
    for(i in 1:n){
    N[i] ~ dnorm(exp.N[i], tau)
    exp.N[i] <- mu.alpha + alpha[site[i]] + beta*TL[i] 
    LogLik[i] <- logdensity.norm(N[i], exp.N[i], tau)
    }
    
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(N = d$Nstd, 
                  TL = d$TotalLength,
                  n = length(d$d15N), 
                  site = as.numeric(site),
                  Nsite = length(unique(site)))


# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         beta = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("beta", "alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")

# Call JAGS from R
mod.N.TLsite = jags(jags.data,
                    inits = NULL,
                    parameters.to.save= jags.params,
                    model.file = "N_v_Length_site.txt",
                    n.chains = nc,
                    n.burnin =nb,
                    n.thin = nt,
                    n.iter = ni,
                    DIC = TRUE)

#check model output
#look for convergence and compare sigma values
mod.N.TLsite
#mcmcplot(mod.N.TLsite)

#calculate total variance
tot_var <- mod.N.TLsite$BUGSoutput$mean$sigma + 
  mod.N.TLsite$BUGSoutput$mean$sigma.alpha

#calculate variance associated with catch site parameter
site_var_ratio <- (mod.N.TLsite$BUGSoutput$mean$sigma.alpha)/tot_var



#plot my results..
MCMCplot(mod.N.TLsite,
         params = c("alpha", "mu.alpha"))
MCMCplot(mod.N.TLsite,
         params = c("beta"))


#plot my results..
ggplot(data = d, aes(x = TotalLength, y = Nstd)) +
  geom_point() +
  geom_abline(intercept = mod.N.TLsite$BUGSoutput$sims.list$mu.alpha[smple], 
              slope = mod.N.TLsite$BUGSoutput$sims.list$beta[smple], 
              color = "gray", alpha = 0.2)+
  geom_abline(slope = mod.N.TLsite$BUGSoutput$mean$beta, 
              intercept = mod.N.TLsite$BUGSoutput$mean$mu.alpha, 
              size  = 1, color = "turquoise4") +
  labs(x = "Total Length (cm)", y = expression(delta~" "^15~N)) +
  theme_classic()


####RANDOM EFFECT OF YEAR####

#Carbon v. Length + (random effect of year)
# Specify model in BUGS language
sink("C_v_Length_year.txt")
cat("
    model {
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of year
    for(j in 1:Nyear){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    beta ~ dnorm(0, 1)
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    
    # Likelihood
    
    for(i in 1:n){
    C[i] ~ dnorm(exp.C[i], tau)
    exp.C[i] <- mu.alpha + alpha[year[i]] + beta*TL[i] 
    LogLik[i] <- logdensity.norm(C[i], exp.C[i], tau)
    }
    
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(C = d$Cstd, 
                  TL = d$TotalLength,
                  n = length(d$d13C), 
                  year = as.numeric(year),
                  Nyear = length(unique(year)))


# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         beta = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("beta", "alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")


# Call JAGS from R
mod.C.TLyr = jags(jags.data,
                  inits = NULL,
                  parameters.to.save= jags.params,
                  model.file = "C_v_Length_year.txt",
                  n.chains = nc,
                  n.burnin =nb,
                  n.thin = nt,
                  n.iter = ni,
                  DIC = TRUE)

#check model output
#look for convergence and compare sigmas
mod.C.TLyr
#mcmcplot(mod.C.TLyr)


#calculate total variance
tot_var <- mod.C.TLyr$BUGSoutput$mean$sigma^2 + 
  mod.C.TLyr$BUGSoutput$mean$sigma.alpha^2

#calculate variance associated with catch site parameter
site_var_ratio <- (mod.C.TLyr$BUGSoutput$mean$sigma.alpha)/tot_var



#plot my results..
MCMCplot(mod.C.TLyr,
         params = c("alpha", "mu.alpha"))
MCMCplot(mod.C.TLyr,
         params = c("beta"))



#plot my results..
ggplot(data = d, aes(x = TotalLength, y = Cstd)) +
  geom_point() +
  geom_abline(intercept = mod.C.TLyr$BUGSoutput$sims.list$mu.alpha[smple], 
              slope = mod.C.TLyr$BUGSoutput$sims.list$beta[smple], 
              color = "gray", alpha = 0.2)+
  geom_abline(slope = mod.C.TLyr$BUGSoutput$mean$beta, 
              intercept = mod.C.TLyr$BUGSoutput$mean$mu.alpha, 
              size  = 1, color = "turquoise4") +
  labs(x = "Total Length (cm)", y = expression(delta~" "^13~C)) +
  theme_classic()

# Now let's do the same model but with Nitrogen.. 

#Nitrogen v. Length + (random effect of year)
sink("N_v_Length_year.txt")
cat("
       model {
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of year
    for(j in 1:Nyear){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    beta ~ dnorm(0, 1)
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    
    # Likelihood
    
    for(i in 1:n){
    N[i] ~ dnorm(exp.N[i], tau)
    exp.N[i] <- mu.alpha + alpha[year[i]] + beta*TL[i] 
    LogLik[i] <- logdensity.norm(N[i], exp.N[i], tau)
    }
    
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(N = d$Nstd, 
                  TL = d$TotalLength,
                  n = length(d$d15N), 
                  year = as.numeric(year),
                  Nyear = length(unique(year)))


# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         beta = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("beta", "alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")

# Call JAGS from R
mod.N.TLyr = jags(jags.data,
                  inits = NULL,
                  parameters.to.save= jags.params,
                  model.file = "N_v_Length_year.txt",
                  n.chains = nc,
                  n.burnin =nb,
                  n.thin = nt,
                  n.iter = ni,
                  DIC = TRUE)


#check model output
#look for convergence and compare sigmas
mod.N.TLyr
#mcmcplot(mod.N.TLyr)


#calculate total variance
tot_var <- mod.N.TLyr$BUGSoutput$mean$sigma + 
  mod.N.TLyr$BUGSoutput$mean$sigma.alpha 

#calculate variance associated with catch site parameter
site_var_ratio <- (mod.N.TLyr$BUGSoutput$mean$sigma.alpha)/tot_var



#plot my results..
MCMCplot(mod.N.TLyr,
         params = c("alpha", "mu.alpha"))
MCMCplot(mod.N.TLyr,
         params = c("beta"))



#convert JAGS output into new format which is easier to plot

#create objects out of output
mn_yr_alpha <- mod.N.TLyr$BUGSoutput$mean$alpha
drws_yr_alpha <- mod.N.TLyr$BUGSoutput$sims.list$alpha

mn_mu_alpha <- mod.N.TLyr$BUGSoutput$mean$mu.alpha
drws_mu_alpha <- mod.N.TLyr$BUGSoutput$sims.list$mu.alpha

mn_beta <- mod.N.TLyr$BUGSoutput$mean$beta
drws_beta <- mod.N.TLyr$BUGSoutput$sims.list$beta

#combine alpha estimates into a single mean estimate for each yr
alpha_muyr <- matrix(NA, nrow = nrow(mn_yr_alpha))
for(i in 1:nrow(mn_yr_alpha)){
  alpha_muyr[i] <- mn_yr_alpha[i] + mn_mu_alpha  
}

#combine for the posterior draws as well
alpha_drwsyr <- matrix(NA, nrow(drws_yr_alpha), 
                         ncol = ncol(drws_yr_alpha))
for(i in 1:ncol(drws_yr_alpha)){
  alpha_drwsyr[,i] <- drws_yr_alpha[,i] + drws_mu_alpha  
}

#assign sample yrs as column names
colnames(alpha_drwsyr) <- c(1:4)

#create a new dataframe for plotting
post_data <- alpha_drwsyr %>% 
  as_tibble(.name_repair = "unique") %>% 
  pivot_longer(cols = 1:4, names_to = "yr", values_to = "alpha_drws")
post_data$alphamn <- rep(alpha_muyr, 39000)
post_data$betadrws <- rep(drws_beta, 4)
post_data$betamn <- rep(mn_beta, nrow(post_data))

#plot my results some more..
colors <- viridis(4)
lsz <- 0.75
smpl <- sample(nrow(post_data), 1000)

d$year <- factor(d$year)

NvLvY <- ggplot() +
  scale_color_manual(values = colors) +
  geom_abline(intercept = post_data$alpha_drws[smpl],
              slope = post_data$betadrws[smpl],
              color = "gray70", alpha = 0.2) +
  geom_point(data = d, aes(x = TotalLength, y = Nstd, color = year)) +
  geom_abline(intercept = post_data$alphamn[1],
              slope = post_data$betamn[1], 
              color = colors[1], size = lsz)+
  geom_abline(intercept = post_data$alphamn[2],
              slope = post_data$betamn[2], 
              color = colors[2], size = lsz)+
  geom_abline(intercept = post_data$alphamn[3],
              slope = post_data$betamn[3], 
              color = colors[3], size = lsz)+
  geom_abline(intercept = post_data$alphamn[4],
              slope = post_data$betamn[4], 
              color = colors[4], size = lsz)+
  labs(x = "Total Length (cm)", y = expression(delta~""^15~N),
       color = "Sample Year") +
  theme_classic()

png(filename="MSFigures/NvTLvYr.png", 
    units="in", 
    width=7, 
    height=6, 
    pointsize=8, 
    res=400)
NvLvY

dev.off()

#Or make everyone's lives easier and just plot mu alpha estimate??
pNvL <- ggplot() +
  geom_abline(intercept = drws_mu_alpha[smple],
              slope = drws_beta[smple],
              color = "gray70", alpha = 0.2) +
  geom_point(data = d, aes(x = TotalLength, y = Nstd)) +
  geom_abline(intercept = mn_mu_alpha,
              slope = mn_beta, 
              size = lsz) +
  labs(x = "Total Length (cm)", y = expression(delta~""^15~N)) +
  theme_classic()

png(filename="MSFigures/NvTL.png", 
    units="in", 
    width=6, 
    height=6, 
    pointsize=8, 
    res=400)
pNvL

dev.off()


####EFFECT OF TL with random effect of YEAR AND CATCHLOCATION####

#Carbon v. Length + (random effect of year/site)
# Specify model in BUGS language
sink("C_v_Length_yrsite.txt")
cat("
    model {
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of year/site
    for(j in 1:Nyrsite){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    beta ~ dnorm(0, 1)
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    
    # Likelihood
    
    for(i in 1:n){
    C[i] ~ dnorm(exp.C[i], tau)
    exp.C[i] <- mu.alpha + alpha[yr_site[i]] + beta*TL[i] 
    LogLik[i] <- logdensity.norm(C[i], exp.C[i], tau)
    }
    
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(C = d$Cstd, 
                  TL = d$TotalLength,
                  n = length(d$d13C), 
                  yr_site = as.numeric(yr_site),
                  Nyrsite = length(unique(yr_site)))


# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         beta = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("beta", "alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")


# Call JAGS from R
mod.C.TLyrsite = jags(jags.data,
                      inits = NULL,
                      parameters.to.save= jags.params,
                      model.file = "C_v_Length_yrsite.txt",
                      n.chains = nc,
                      n.burnin =nb,
                      n.thin = nt,
                      n.iter = ni,
                      DIC = TRUE)


#check model output
#look for convergence and compare sigmas
mod.C.TLyrsite
#mcmcplot(mod.C.TLyrsite)


#calculate total variance
tot_var <- mod.C.TLyrsite$BUGSoutput$mean$sigma + 
  mod.C.TLyrsite$BUGSoutput$mean$sigma.alpha

#calculate variance associated with catch site parameter
site_var_ratio <- (mod.C.TLyrsite$BUGSoutput$mean$sigma.alpha)/tot_var



#plot my results..
MCMCplot(mod.C.TLyrsite,
         params = c("alpha", "mu.alpha"))
MCMCplot(mod.C.TLyrsite,
         params = c("beta"))


#plot my results..
ggplot(data = d, aes(x = TotalLength, y = Cstd)) +
  geom_point() +
  geom_abline(intercept = mod.C.TLyrsite$BUGSoutput$sims.list$mu.alpha[smple], 
              slope = mod.C.TLyrsite$BUGSoutput$sims.list$beta[smple], 
              color = "gray", alpha = 0.2)+
  geom_abline(slope = mod.C.TLyrsite$BUGSoutput$mean$beta, 
              intercept = mod.C.TLyrsite$BUGSoutput$mean$mu.alpha, 
              size  = 1, color = "turquoise4") +
  labs(x = "Total Length (cm)", y = expression(delta~" "^13~C)) +
  theme_classic()


# Now let's do the same model for N! 
#N v. Length + (random effect of year and catchsite)
# Specify model in BUGS language
sink("N_v_Length_yrsite.txt")
cat("
    model {
    
    # Priors
    mu.alpha ~ dnorm(0, 1) #mean hyperparameter for random intercepts
    sigma.alpha ~ dunif(0, 5) #SD hyperparameter for random intercepts
    tau.alpha <- 1/(sigma.alpha*sigma.alpha) #convert to tau
    
    
    # Random effect of year/site
    for(j in 1:Nyrsite){
    alpha[j] ~ dnorm(0, tau.alpha) #random intercepts
    }
    
    beta ~ dnorm(0, 1)
    
    sigma ~ dexp(1) #residual standard deviation
    tau <- 1/(sigma*sigma) #convert to tau
    
    
    # Likelihood
    
    for(i in 1:n){
    N[i] ~ dnorm(exp.N[i], tau)
    exp.N[i] <- mu.alpha + alpha[yr_site[i]] + beta*TL[i] 
    LogLik[i] <- logdensity.norm(N[i], exp.N[i], tau)
    }
    
    
    }
    ", fill = TRUE)
sink()


# Bundle data
jags.data <- list(N = d$Nstd, 
                  TL = d$TotalLength,
                  n = length(d$d15N), 
                  yr_site = as.numeric(yr_site),
                  Nyrsite = length(unique(yr_site)))


# Initial values
inits <- function() list(alpha = rnorm(1, 0, 0.5), 
                         mu.alpha = rnorm(1, 0, 0.5),
                         beta = rnorm(1, 0, 0.5),
                         sigma = runif(1, 0, 1),
                         sigma.alpha = runif(1, 0, 1))

# Parameters you want to keep track of:
jags.params <- c("beta", "alpha", "mu.alpha", "sigma", "sigma.alpha", "LogLik")


# Call JAGS from R
mod.N.TLyrsite = jags(jags.data,
                      inits = NULL,
                      parameters.to.save= jags.params,
                      model.file = "N_v_Length_yrsite.txt",
                      n.chains = nc,
                      n.burnin =nb,
                      n.thin = nt,
                      n.iter = ni,
                      DIC = TRUE)


#check model output
#look for convergence and compare sigmas
mod.N.TLyrsite
#mcmcplot(mod.N.TLyrsite)

#calculate total variance
tot_var <- mod.N.TLyrsite$BUGSoutput$mean$sigma + 
  mod.N.TLyrsite$BUGSoutput$mean$sigma.alpha

#calculate variance associated with catch site parameter
site_var_ratio <- (mod.N.TLyrsite$BUGSoutput$mean$sigma.alpha)/tot_var



#plot my results..
MCMCplot(mod.N.TLyrsite,
         params = c("alpha", "mu.alpha"))
MCMCplot(mod.N.TLyrsite,
         params = c("beta"))


#plot my results..
ggplot(data = d, aes(x = TotalLength, y = Nstd)) +
  geom_point() +
  geom_abline(intercept = mod.N.TLyrsite$BUGSoutput$sims.list$mu.alpha[smple], 
              slope = mod.N.TLyrsite$BUGSoutput$sims.list$beta[smple], 
              color = "gray", alpha = 0.2)+
  geom_abline(slope = mod.N.TLyrsite$BUGSoutput$mean$beta, 
              intercept = mod.N.TLyrsite$BUGSoutput$mean$mu.alpha, 
              size  = 1, color = "turquoise4") +
  labs(x = "Total Length (cm)", y = expression(delta~" "^15~N)) +
  theme_classic()


####Model Comparison/Selection####

#CARBON

CTL <- loo(mod.C.TL$BUGSoutput$sims.list$LogLik)
CSite <- loo(mod.C.site$BUGSoutput$sims.list$LogLik)
Cyr <- loo(mod.C.year$BUGSoutput$sims.list$LogLik)
Cyrsite <- loo(mod.C.yrsite$BUGSoutput$sims.list$LogLik)
CTLsite <- loo(mod.C.TLsite$BUGSoutput$sims.list$LogLik)
CTLyr <- loo(mod.C.TLyr$BUGSoutput$sims.list$LogLik)
CTLyrsite <- loo(mod.C.TLyrsite$BUGSoutput$sims.list$LogLik)


loo_compare(CTL, CSite, Cyr, Cyrsite, CTLsite, CTLyr, CTLyrsite)
#looks like carbon v. total length is the best model?

#let's just see what waic has to say
CTL.w <- waic(mod.C.TL$BUGSoutput$sims.list$LogLik)
CSite.w <- waic(mod.C.site$BUGSoutput$sims.list$LogLik)
Cyr.w <- waic(mod.C.year$BUGSoutput$sims.list$LogLik)
Cyrsite.w <- waic(mod.C.yrsite$BUGSoutput$sims.list$LogLik)
CTLsite.w <- waic(mod.C.TLsite$BUGSoutput$sims.list$LogLik)
CTLyr.w <- waic(mod.C.TLyr$BUGSoutput$sims.list$LogLik)
CTLyrsite.w <- waic(mod.C.TLyrsite$BUGSoutput$sims.list$LogLik)

loo_compare(CTL.w, CSite.w, Cyr.w, Cyrsite.w, CTLsite.w, CTLyr.w, CTLyrsite.w)
#this says CTLyrsite...

#both methods indicate TL, site, and yrsite do a similarly good job at predicting the relationship

#NITROGEN
NTL <- loo(mod.N.TL$BUGSoutput$sims.list$LogLik)
NSite <- loo(mod.N.site$BUGSoutput$sims.list$LogLik)
Nyr <- loo(mod.N.year$BUGSoutput$sims.list$LogLik)
Nyrsite <- loo(mod.N.yrsite$BUGSoutput$sims.list$LogLik)
NTLsite <- loo(mod.N.TLsite$BUGSoutput$sims.list$LogLik)
NTLyr <- loo(mod.N.TLyr$BUGSoutput$sims.list$LogLik)
NTLyrsite <- loo(mod.N.TLyrsite$BUGSoutput$sims.list$LogLik)


loo_compare(NTL, NSite, Nyr, Nyrsite, NTLsite, NTLyr, NTLyrsite)
#looks like N v. year is the best model?

#let's just see what waic has to say
NTL.w <- waic(mod.N.TL$BUGSoutput$sims.list$LogLik)
NSite.w <- waic(mod.N.site$BUGSoutput$sims.list$LogLik)
Nyr.w <- waic(mod.N.year$BUGSoutput$sims.list$LogLik)
Nyrsite.w <- waic(mod.N.yrsite$BUGSoutput$sims.list$LogLik)
NTLsite.w <- waic(mod.N.TLsite$BUGSoutput$sims.list$LogLik)
NTLyr.w <- waic(mod.N.TLyr$BUGSoutput$sims.list$LogLik)
NTLyrsite.w <- waic(mod.N.TLyrsite$BUGSoutput$sims.list$LogLik)

loo_compare(NTL.w, NSite.w, Nyr.w, Nyrsite.w, NTLsite.w, NTLyr.w, NTLyrsite.w)


###ulam models for NTLyr and CTLsite
library(rethinking)

d2 <- select(d, Nstd, year, TotalLength)
d2$year <- as.numeric(as.factor(d2$year))


m.NTLyr <- ulam(
  alist(
        Nstd ~ dnorm(expN, sigma),
        expN <- a[year] + beta*TotalLength,
        a[year] ~ dnorm(a_bar, sigma_bar),
        a_bar ~ dnorm(0, 1),
        sigma_bar ~ dexp(1),
        sigma ~ dexp(1),
        beta ~ dnorm(0,1)
      ), data=d2 , chains=4 ,  log_lik=TRUE )

d2a <- select(d, Nstd, TotalLength, yr_site)
d2a$yr_site <- as.numeric(as.factor(d2a$yr_site))

m.NTLyrsite <- ulam(
  alist(
    Nstd ~ dnorm(expN, sigma),
    expN <- a[yr_site] + beta*TotalLength,
    a[yr_site] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1),
    beta ~ dnorm(0,1)
  ), data=d2a , chains=4 ,  log_lik=TRUE )


d2b <- select(d, Nstd, TotalLength, sample_loc)
d2b$site <- as.numeric(as.factor(d2b$sample_loc))

m.TLsite <- ulam(
  alist(
    Nstd ~ dnorm(expN, sigma),
    expN <- a[site] + beta*TotalLength,
    a[site] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1),
    beta ~ dnorm(0,1)
  ), data=d2b , chains=4 ,  log_lik=TRUE )


m.site <- ulam(
  alist(
    Nstd ~ dnorm(expN, sigma),
    expN <- a[site],
    a[site] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1)
  ), data=d2b , chains=4 ,  log_lik=TRUE 
)

m.yr <- ulam(
  alist(
    Nstd ~ dnorm(expN, sigma),
    expN <- a[year],
    a[year] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1)
  ), data=d2 , chains=4 ,  log_lik=TRUE 
)

m.yrsite <- ulam(
  alist(
    Nstd ~ dnorm(expN, sigma),
    expN <- a[yr_site],
    a[yr_site] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1)
  ), data=d2a , chains=4 ,  log_lik=TRUE 
)

m.TL <- ulam(
  alist(
    Nstd ~ dnorm(expN, sigma),
    expN <- a + beta*TotalLength,
    a ~ dnorm(0, 1),
    sigma ~ dexp(1),
    beta ~ dnorm(0,1)
  ), data=d2 , chains=4 ,  log_lik=TRUE )

compare(m.NTLyr, m.NTLyrsite, m.TLsite, m.site, m.yr, m.yrsite, m.TL)


#Carbon
d3 <- select(d, Cstd, year, TotalLength)
d3$year <- as.numeric(as.factor(d3$year))


m.CTLyr <- ulam(
  alist(
    Cstd ~ dnorm(expC, sigma),
    expC <- a[year] + beta*TotalLength,
    a[year] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1),
    beta ~ dnorm(0,1)
  ), data=d3 , chains=4 ,  log_lik=TRUE )

d3a <- select(d, Cstd, TotalLength, yr_site)
d3a$yr_site <- as.numeric(as.factor(d3a$yr_site))

m.CTLyrsite <- ulam(
  alist(
    Cstd ~ dnorm(expC, sigma),
    expC <- a[yr_site] + beta*TotalLength,
    a[yr_site] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1),
    beta ~ dnorm(0,1)
  ), data=d3a , chains=4 ,  log_lik=TRUE )


d3b <- select(d, Cstd, TotalLength, sample_loc)
d3b$site <- as.numeric(as.factor(d3b$sample_loc))

m.TLsite <- ulam(
  alist(
    Cstd ~ dnorm(expC, sigma),
    expC <- a[site] + beta*TotalLength,
    a[site] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1),
    beta ~ dnorm(0,1)
  ), data=d3b , chains=4 ,  log_lik=TRUE )


m.Csite <- ulam(
  alist(
    Cstd ~ dnorm(expC, sigma),
    expC <- a[site],
    a[site] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1)
  ), data=d3b , chains=4 ,  log_lik=TRUE 
)

m.Cyr <- ulam(
  alist(
    Cstd ~ dnorm(expC, sigma),
    expC <- a[year],
    a[year] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1)
  ), data=d3 , chains=4 ,  log_lik=TRUE 
)

m.Cyrsite <- ulam(
  alist(
    Cstd ~ dnorm(expC, sigma),
    expC <- a[yr_site],
    a[yr_site] ~ dnorm(a_bar, sigma_bar),
    a_bar ~ dnorm(0, 1),
    sigma_bar ~ dexp(1),
    sigma ~ dexp(1)
  ), data=d3a , chains=4 ,  log_lik=TRUE 
)

m.CTL <- ulam(
  alist(
    Cstd ~ dnorm(expC, sigma),
    expC <- a + beta*TotalLength,
    a ~ dnorm(0, 1),
    sigma ~ dexp(1),
    beta ~ dnorm(0,1)
  ), data=d3 , chains=4 ,  log_lik=TRUE )

compare(m.CTLyr, m.CTLyrsite, m.TLsite, m.Csite, m.Cyr, m.Cyrsite, m.CTL)
