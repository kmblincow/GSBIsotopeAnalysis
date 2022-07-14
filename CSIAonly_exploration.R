#Kayla Blincow

#exploring the CSIA samples 

rm(list = ls())


library(tidyverse)
library(R2jags)



#load the data
gsb <- read.csv("CSIA_length.csv")


gsb <- filter(gsb, !Phe > 15 & !is.na(TP_Fcerror))

#look for relationship between TP and TL

sink("TP_v_TL.txt")
cat("
    model {
    
    # Priors
    
    alpha ~ dunif(0, 50)
    beta ~ dunif(0,5)
    
    sigma ~ dunif(0,10)
    tau <- 1/(sigma*sigma)
    
    # Likelihood
    
    for(i in 1:n){
    exp.TP[i] <- alpha + beta * TL[i]
    TP[i] ~ dnorm(exp.TP[i], tau)    
    }
    
    }
    ", fill = TRUE)
sink()

# Bundle data
jags.data <- list(TP = gsb$TP_Fcerror, n = length(gsb$BulkN), 
                  TL = gsb$TotalLength)

# Initial values
inits <- function() list(alpha = runif(1, 0, 10), sigma = runif(1, 0, 5), beta = runif(1, -10, 10))

# Parameters you want to keep track of:
jags.params <- c("alpha", "beta", "sigma")

# MCMC settings
ni <- 350000 #iterations
nt <- 25 #thinning by
nb <- 25000 #burn in
nc <- 3 #chains

# Call JAGS from R
mod = jags(jags.data,
           inits = NULL,
           parameters.to.save= jags.params,
           model.file= "TP_v_TL.txt",
           n.chains = nc,
           n.burnin =nb,
           n.thin = nt,
           n.iter = ni,
           DIC = TRUE)

#check out the model
mod


#generate 500 values to create index to randomly sample posterior
smple <- sample(1:39000, 1000)

#plot my results..
pTP <- ggplot(data = gsb, aes(x = TotalLength, y = TP_Fcerror)) +
  geom_abline(intercept = mod$BUGSoutput$sims.list$alpha[smple],
              slope = mod$BUGSoutput$sims.list$beta[smple],
              color = "gray", alpha = 0.2)+
  geom_abline(slope = mod$BUGSoutput$mean$beta,
              intercept = mod$BUGSoutput$mean$alpha,
              size  = 1, color = "black") +
  geom_point() +
  labs(x = "Total Length (cm)", y = "CSIA-AA Trophic Position Estimate") +
  theme_classic()

png(filename="C:/Users/kmbli/OneDrive - UC San Diego/PhDzNuts/Giant Sea Bass/Isotopes/IsotopesAnalysis/MSFigures/TPvLength.png", 
    units="in", 
    width=6, 
    height=6, 
    pointsize=8, 
    res=400)
pTP

dev.off()


#plots of the the AA values
library(patchwork)
d <- read.csv("AAvalues.csv", header = T)
d <- filter(d, !Phe > 15)

d2 <- pivot_longer(d, cols = Ala:Phe, names_to = "AA", values_to = "Value")
d2$AA <- factor(d2$AA, labels = c("Alanine", "Glutamic Acid", "Glycine", "Phenylalanine"))

aa <- ggplot(d2) +
  geom_point(aes(x = TotalLength, y = Value, shape = AA), size = 2) +
  labs(x = "Total Length (cm)", y = expression(delta~" "^15~N), 
       shape = "Amino Acid") +
  theme_classic()

png(filename="C:/Users/kmbli/OneDrive - UC San Diego/PhDzNuts/Giant Sea Bass/Isotopes/IsotopesAnalysis/MSFigures/AAPlot.png", 
    units="in", 
    width=7, 
    height=6, 
    pointsize=8, 
    res=400)
aa

dev.off()


#Outlier removal supplement
library(ggforce)
d <- read.csv("AAvalues.csv", header = T)
d$color <- "Black"

d[d$ArtID == "EN 170610-5" | d$ArtID == "EN 170701-1" | d$ArtID == "EN 170713-1",]$color <- "Red"

d2 <- pivot_longer(d, cols = Ala:Phe, names_to = "AA", values_to = "Value")
d2$AA <- factor(d2$AA, labels = c("Alanine", "Glutamic Acid", "Glycine", "Phenylalanine"))


outlier <- ggplot(d2) +
  geom_point(aes(x = TotalLength, y = Value, shape = AA, color = color), 
             size = 2) +
  scale_color_manual(values = c("black", "red")) +
  geom_ellipse(aes(x0 = 75, y0 = 16.1, a = 30, b = 1, angle = 0),
               color = "red", size = 1) +
  labs(x = "Total Length (cm)", y = expression(delta~" "^15~N), 
       shape = "Amino Acid") +
  guides(color = FALSE) +
  theme_classic()

png(filename="C:/Users/kmbli/OneDrive - UC San Diego/PhDzNuts/Giant Sea Bass/Isotopes/IsotopesAnalysis/MSFigures/OutlierRemoval.png", 
    units="in", 
    width=7, 
    height=6, 
    pointsize=8, 
    res=400)
outlier

dev.off()



#Other stuff

all <- read.csv("FinalGSBBulk.csv")
ggplot(all, aes(x = TotalLength, y = d15N)) +
  geom_point() +
  stat_smooth(method = "lm") +
  facet_wrap(vars(sample_loc))

ggplot(gsb, aes(x = TotalLength, y = d13C)) +
  geom_point() +
  stat_smooth(method = "lm")


ggplot(filter(all, sample_loc == "GuerreroNegro"), 
       aes(x = TotalLength, y = d15N)) +
  stat_smooth(method = "lm") +
  geom_point()

