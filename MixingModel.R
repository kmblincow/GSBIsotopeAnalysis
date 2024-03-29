#Kayla Blincow
#Update: 3/2/2021

#The goal of this script is to run a Bayesian mixing model to look at relative
#contributions of different primary production sources to GSB diets using MixSIAR.

#Basing this off the alligator example from Stock PeerJ publication

#Update 3/9/2021: after talking to Brice decided adding a random effect of individual
#would confounded with individual effects associated with length fixed effect

#clear my workspace 
rm(list = ls())

#load packages
library(tidyverse)
library(MixSIAR)
library(R2jags)



####Run CSIA derived N TDF mixing model####

# load mixture (consumer) data

mix <- load_mix_data(filename="FinalGSBBulk.csv",
                     iso_names=c("d13C", "d15N"),
                     factors = NULL,
                     fac_random=NULL,
                     fac_nested=NULL,
                     cont_effects="TotalLength")

# mix <- load_mix_data(filename="CSIA_length.csv",
#                      iso_names=c("d13C", "d15N"),
#                      factors = NULL,
#                      fac_random=NULL,
#                      fac_nested=NULL,
#                      cont_effects="TotalLength")


#with fixed effect of sample size
source <- load_source_data(filename= "MixingModel/IsotopeMixingModel_GSB/IsotopeMixingModel_GSB/PPSources_means.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix=mix)

# load TEF data
#TEF data with calculated slope value for N and standard literature values for C
discr <- load_discr_data(filename = "MixingModel/IsotopeMixingModel_GSB/IsotopeMixingModel_GSB/gsb_TEF_stCrnge.csv", mix = mix)


#isospace plot
plot_data(filename="isospace_plot",
          plot_save_pdf=FALSE,
          plot_save_png=TRUE,
          mix, source, discr)


# Define model structure and write JAGS model file
model_filename <- paste0("MixSIAR_model_length_stdCTEF.txt")
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model
jags.mod <- run_model(run="short", mix, source, discr, model_filename, 
                      alpha.prior=1, resid_err, process_err)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.mod, mix, source,
            output_options = list(
              summary_save = TRUE,
              summary_name = "ModelSummary/summary_statistics_HighSS",
              sup_post = TRUE,
              plot_post_save_pdf = TRUE,
              plot_post_name = "ModelSummary/Posteriors/posterior_density_HighSS",
              sup_pairs = TRUE,
              plot_pairs_save_pdf = TRUE, 
              plot_pairs_name = "ModelSummary/pairs_plot_HighSS",
              sup_xy = TRUE,
              plot_xy_save_pdf = TRUE,
              plot_xy_name = "ModelSummary/traceplot_HighSS",
              gelman = TRUE,
              heidel = FALSE,
              geweke = FALSE,
              diag_save = TRUE,
              diag_name = "ModelSummary/diagnostics_HighSS",
              indiv_effect = FALSE,
              plot_post_save_png = FALSE,
              plot_pairs_save_png = FALSE,
              plot_xy_save_png = FALSE)
)
graphics.off()

####plot the results in terms of proportion of PP source by length####
# Plot proportions vs. length 
R2jags::attach.jags(jags.mod)
n.sources <- source$n.sources
source_names <- source$source_names

get_high <- function(x){return(quantile(x,.975))}
get_low <- function(x){return(quantile(x,.025))}

n.plot = 200

chain.len = dim(p.global)[1]

Cont1.plot <- seq(from = round(min(mix$CE[[1]]), 1), 
                  to = round(max(mix$CE[[1]]), 1), 
                  length.out = n.plot)

ilr.plot <- array(NA, dim = c(n.plot, n.sources-1, chain.len))
ilr.median <- array(NA, dim = c(n.plot, n.sources-1))
ilr.low <- array(NA, dim = c(n.plot, n.sources-1))
ilr.high <- array(NA, dim = c(n.plot, n.sources-1))

for(src in 1:n.sources-1){
  for(i in 1:n.plot){
    ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont1[,src]*Cont1.plot[i]
    ilr.low[i,src] <- get_low(ilr.plot[i,src,])
    ilr.median[i,src] <- median(ilr.plot[i,src,])
    ilr.high[i,src] <- get_high(ilr.plot[i,src,])
  }
}

# Transform regression lines from ILR-space to p-space
e <- matrix(rep(0,n.sources*(n.sources-1)),
            nrow = n.sources, ncol = (n.sources-1))

for(i in 1:(n.sources - 1)){
  e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
  e[,i] <- e[,i]/sum(e[,i])
}

#create a bunch of empty arrays that we will fill with a for loop in a sec
cross.med <- array(data = NA, dim = c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation

tmp.p.med <- array(data = NA, dim = c(n.plot, n.sources))              # dummy variable for inverse ILR calculation

p.median <- array(data = NA, dim = c(n.plot, n.sources))

cross.low <- array(data = NA, dim = c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation

tmp.p.low <- array(data = NA, dim = c(n.plot, n.sources))              # dummy variable for inverse ILR calculation

p.low <- array(data = NA, dim = c(n.plot, n.sources))

cross.high <- array(data = NA, dim = c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation

tmp.p.high <- array(data = NA, dim = c(n.plot, n.sources))              # dummy variable for inverse ILR calculation

p.high <- array(data = NA, dim = c(n.plot, n.sources))

eps.low <- rep(NA, n.plot)
eps.med <- rep(NA, n.plot)
eps.high <- rep(NA, n.plot)    

for(i in 1:n.plot){
  for(j in 1:(n.sources-1)){
    cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
    cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
    cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
  }
  for(src in 1:n.sources){
    tmp.p.med[i,src] <- prod(cross.med[i,src,]);
    tmp.p.low[i,src] <- prod(cross.low[i,src,]);
    tmp.p.high[i,src] <- prod(cross.high[i,src,]);
  }
  for(src in 1:n.sources){
    p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
    p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
    p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
  }
}

colnames(p.median) <- source_names

Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), 
                 reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
colnames(df) <- c("source","median","x","low","high")

# Plot of Diet vs. Cont effect

p <- print(ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = median)) +
             ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high, 
                                               group = source), 
                                  alpha=0.35, fill = "gray60") +
             ggplot2::geom_line(ggplot2::aes(x = x, y = median,
                                             group = source, 
                                             linetype = source), 
                                color = "black", size=1.5) +
             ggplot2::ylab("Diet Proportion") +
             ggplot2::xlab("Total Length (cm)") +
             ggplot2::scale_linetype_discrete(labels = c("Macroalgae",
                                                         "Phytoplankton")) +
             ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
             ggplot2::theme_bw() +
             ggplot2::theme(panel.border = ggplot2::element_blank(), 
                            panel.grid.major = ggplot2::element_blank(), 
                            panel.grid.minor = ggplot2::element_blank(), 
                            panel.background = ggplot2::element_blank(), 
                            axis.line = ggplot2::element_line(colour = "black"), 
                            axis.title = ggplot2::element_text(size = 16), 
                            axis.text = ggplot2::element_text(size = 14), 
                            legend.text = ggplot2::element_text(size = 14), 
                            legend.position = c(.7,.5), 
                            legend.justification = c(0,1), 
                            legend.title = ggplot2::element_blank()))


png("Length_CN_CSIATDF.png", height=7, width=7, units='in', res=500)

p 

dev.off()

#for publication
tiff("MSFigures/FinalFigs/M14123_Fig4.tif",
     width = 2800, 
     height = 2600,
     res = 300)

p

dev.off()
