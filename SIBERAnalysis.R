#Kayla Blincow
#10/29/2020

#The purpose of this script is generate a C/N biplot with standard ellipse areas 
#with each GSB individual grouped based on the age classifications

#clear my workspace
rm(list = ls())

#load my libraries
library(tidyverse)
library(SIBER)
library(viridis)
library(coda)

#load my data
d <- read.csv("FinalGSBBulk.csv", header = T)


#create a colun that's age group
d$group <- d$sclass


#convert to a dataframe that SIBER will accept
d_SIBER <- select(d, d13C, d15N, group) %>% 
  mutate(community = 1)
names(d_SIBER) <- c("iso1", "iso2", "group", "community")

#create a SIBER object
dsib <- createSiberObject(d_SIBER)

#do the plots?
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

palette(viridis(3))


par(mfrow=c(1,1))

png(filename="MSFigures/SIBER_Ellipses.png", 
    units="in", 
    width=4, 
    height=4, 
    pointsize=8, 
    res=400)

plotSiberObject(dsib, ax.pad = 2,
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C),
                ylab = expression({delta}^15*N),
                points.order = 20
)
legend("topright", legend = c("Mature", "Transition", "Immature"), 
       fill = palette(viridis(3)))


plotGroupEllipses(dsib, n = 100, p = NULL,
                  lty = 1, lwd = 2, small.sample = T)

plotGroupEllipses(dsib, n = 100, p = 0.95, small.sample = T,
                  lty = 2, lwd = 1)

dev.off()

#calculate summary statistics
group.ML <- groupMetricsML(dsib)
print(group.ML)


#####Fit Bayesian Models to the Data####
# options for running jags
parms <- list()
parms$n.iter <- 50000   # number of iterations to run the model for
parms$n.burnin <- 25000 # discard the first set of values
parms$n.thin <- 5      # thin the posterior by this many
parms$n.chains <- 3        # run this many chains

# set save.output = TRUE
parms$save.output = TRUE
# you might want to change the directory to your local directory or a 
# sub folder in your current working directory. I have to set it to a 
# temporary directory that R creates and can use for the purposes of this 
# generic vignette that has to run on any computer as the package is 
# built and installed.
parms$save.dir = tempdir()

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(dsib, parms, priors)

#test for convergence
# get a list of all the files in the save directory
all.files <- dir(parms$save.dir, full.names = TRUE)
# find which ones are jags model files
model.files <- all.files[grep("jags_output", all.files)]
# test convergence for the first one
do.this <- 1
load(model.files[do.this])
gelman.diag(output, multivariate = FALSE)
gelman.plot(output, auto.layout = FALSE)
#looks like we've converged!


#plot posteriors for SEAb estimates
# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

#reorganize columns for plotting
col.order <- c(2,1,3)
SEA.B <- SEA.B[,col.order]
group.ML <- group.ML[,col.order]

png(filename="MSFigures/SEA_estimates.png", 
    units="in", 
    width=4, 
    height=4, 
    pointsize=8, 
    res=400)

p <- siberDensityPlot(SEA.B, xticklabels = c("Immature", "Transition", "Mature"), 
                 xlab = c("Age Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 ylims = c(0.25, 6)
)


dev.off()

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)









#look at overlap of ellipses
ellipse1 <- "1.Large"

ellipse2 <- "1.Medium"

ellipse3 <- "1.Small"

bayes95.overlap <- bayesianOverlap(ellipse1, ellipse2, ellipses.posterior,
                                   draws = 100, p.interval = 0.95, n = 100)
# a histogram of the overlap
hist(bayes95.overlap[,3], 10)
# and as above, you can express this a proportion of the non-overlapping area of 
# the two ellipses, would be
bayes.prop.95.over <- (bayes95.overlap[,3] / (bayes95.overlap[,2] + 
                                                bayes95.overlap[,1] -
                                                bayes95.overlap[,3])
)
hist(bayes.prop.95.over, 10)




#####Plot posterior ellipses####
# how many of the posterior draws do you want?
n.posts <- 20

# for a standard ellipse use
p.ell <- pchisq(1,2)

# a list to store the results
all_ellipses <- list()
# loop over groups
for (i in 1:length(ellipses.posterior)){
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  for ( j in 1:n.posts){
    
    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)
    
    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]
    
    # ellipse points
    
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}
ellipse_df <- bind_rows(all_ellipses, .id = "id")
# now we need the group and community names
# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]
# split them and conver to a matrix, NB byrow = T
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)
ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]
ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)

#plot the data
ggplot(d_SIBER) +
  geom_point(aes(iso1, iso2, color = factor(group)), size = 2)+
  geom_polygon(data = ellipse_df,
               mapping = aes(iso1, iso2,
                             color = group),
               alpha = 0.2) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=15)) +
  theme_classic()


