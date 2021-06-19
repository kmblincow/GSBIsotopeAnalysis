#Kayla Blincow
#3/22/2021

#The purpose of this script is calculate TDF values based on the CSIA TP estimates

#Update: 3/30/2021
#Going to calculate the TDF by calculating the slope of the TP vs Bulk N linear
#relationship
#data needed, TP estimates, and bulk N values

#UPDATE: 4/8/2021
#Should not propagate error of TP estimates 

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

#load the data
gsb <- read.csv("CSIA_length.csv")

#remove three outliers (really high Phe values and really low TP values)
gsb <- gsb[!gsb$Phe > 15 & !is.na(gsb$TP_Fcerror),]

#plot the relationship
plot(gsb$TP_Fcerror, gsb$BulkN)

#fit a linear model in JAGS

#Bulk N v. TP#
# Specify model in BUGS language
sink("N_v_TP.txt")
cat("
    model {
    
    # Priors
    
    alpha ~ dunif(0, 50)
    beta ~ dunif(0,5)
    
    sigma ~ dunif(0,10)
    tau <- 1/(sigma*sigma)
    
    # Likelihood
    
    for(i in 1:n){
    exp.N[i] <- alpha + beta * TP[i]
    N[i] ~ dnorm(exp.N[i], tau)    
    }
    
    }
    ", fill = TRUE)
sink()

# Bundle data
jags.data <- list(N = gsb$BulkN, n = length(gsb$BulkN), TP = gsb$TP_Fcerror)

# Initial values
inits <- function() list(alpha = runif(1, 0, 10), sigma = runif(1, 0, 5), beta = runif(1, -10, 10))

# Parameters you want to keep track of:
jags.params <- c("alpha", "beta", "sigma")

# MCMC settings
ni <- 200000 #iterations
nt <- 10 #thinning by
nb <- 25000 #burn in
nc <- 3 #chains

# Call JAGS from R
mod = jags(jags.data,
           inits = NULL,
           parameters.to.save= jags.params,
           model.file= "N_v_TP.txt",
           n.chains = nc,
           n.burnin =nb,
           n.thin = nt,
           n.iter = ni,
           DIC = TRUE)

#check out the model
mod

MCMCplot(mod, params = "beta")

#generate 500 values to create index to randomly sample posterior
smple <- sample(1:52500, 1000)

#plot my results..
ggplot(data = gsb, aes(x = TP_Fcerror, y = BulkN)) +
  geom_point() +
  geom_abline(intercept = mod$BUGSoutput$sims.list$alpha[smple], 
              slope = mod$BUGSoutput$sims.list$beta[smple], 
              color = "gray", alpha = 0.2)+
  geom_abline(slope = mod$BUGSoutput$mean$beta, 
              intercept = mod$BUGSoutput$mean$alpha, 
              size  = 1, color = "turquoise4") +
  labs(x = "Trophic Position", y = expression(delta~" "^15~N)) +
  theme_classic()

#calculate the appropriate N TDF for mixing model (slope * mean(TP))
NTDF <- mod$BUGSoutput$mean$beta*(mean(gsb$TP_Fcerror)) - 1
#1.57 is the mean TDF from primary producers to GSBs 
sd_NTDF <- mod$BUGSoutput$sd$beta
#0.39

#calculate N TDF from standard 3.4 value * mean(TP)
NTDF_st <- 3.4*(mean(gsb$TP_Fcerror)) - 1 #10.538
#sd = 1



#caclulate the C TDF value from reasonable value based onMichener and Kaufman review (0.9 +/- 0.5)

CTDF_std <- 0.9*(mean(gsb$TP_Fcerror)) - 1
#2.054

#sd = 0.5

####OLD WAY####
#data needed: CSIA data, bulk N/C values for primary producers

#clear my workspace
rm(list = ls())


#load the data
gsb <- read.csv("CSIA_length.csv")
pp <- read.csv("PPsources.csv")

#write the equation that doesn't account for different sources
TDF <- function(N15GSB, N15PP, TPGSB){
  tdf = (N15GSB - N15PP)/(TPGSB - 1)
  print(tdf)
}

#calculate the TDF including mean and sd (for mean of pp options)
tdf1 <- tdf1a <- TDF(N15GSB = gsb$BulkN, N15PP = mean(pp$Meand15N), TPGSB = gsb$TP_Fcerror)
mean(tdf1, na.rm = T) #4.476189
sd(tdf1, na.rm = T) #0.746189

#values input into gsb_TEF_mean.csv



#calculate the TDF including mean and sd (for each pp option)
tdf1a <- TDF(N15GSB = gsb$BulkN, N15PP = pp$Meand15N[1], TPGSB = gsb$TP_Fcerror)
mean(tdf1a, na.rm = T) #3.857867 (for KLP)
sd(tdf1a, na.rm = T) #0.6399658 (for KLP)

tdf1b <- TDF(N15GSB = gsb$BulkN, N15PP = pp$Meand15N[2], TPGSB = gsb$TP_Fcerror)
mean(tdf1b, na.rm = T) #5.086456 (for PHYTO)
sd(tdf1b, na.rm = T) #0.8541674 (for PHYTO)

#values input into gsb_TEF_sep.csv



#write equation that does account for different sources
TDF2 <- function(C13GSB, C13PP1, C13PP2, N15GSB, N15PP1, N15PP2, TPGSB){
  alpha = (C13GSB - C13PP2)/(C13PP1 - C13PP2)
  tdf = 1 + (N15GSB - (N15PP1 * alpha + N15PP2 * (1-alpha)))/TPGSB
  print(tdf)
}

#calculate the TDF including mean and sd
tdf2 <- TDF2(C13GSB = gsb$BulkC, C13PP1 = pp$Meand13C[1], C13PP2 = pp$Meand13C[2],
             N15GSB = gsb$BulkN, N15PP1 = pp$Meand15N[1], N15PP2 = pp$Meand15N[2],
             TPGSB = gsb$TP_Fcerror)
mean(tdf2, na.rm = T) #3.61754
sd(tdf2, na.rm = T) #0.4304933
#values input into gsb_TEF_mix.csv