#Kayla Blincow
#3/9/2021

#Base Plots for GSB Isotope Analysis


#The purpose of this script is to generate all the plots/analysis for the basis
#of the GSB isotope chapter.
#Much of this code exists elsewhere, but it's all over the place. This is your
#one stop shop for the plots with all of the bulk data.

#clear my workspace
rm(list = ls())

#load packages
library(tidyverse)

#load ze data!
d <- read.csv("FinalGSBBulk.csv", header = T)


#####Isotope Bi-Plot (with catchsite and length dimension)####
biplot <- ggplot(d, aes(x = d13C, y = d15N, size = TotalLength)) +
  geom_point() +
  labs(x = expression(delta~""^13~C),
       y = expression(delta~""^15~N),
       shape = "Capture Location",
       size = "Total Length (cm)") +
  theme_classic() 
  
png(filename="MSFigures/Biplot.png", 
    units="in", 
    width=7, 
    height=6, 
    pointsize=8, 
    res=400)

biplot

dev.off()

####Bayesian Linear Models to test c/N v length, catchsite, and length + catchsite####
#This information is not located in the "JAGS_linearmodel.R" script

#plots of catch location information to ask question.
ggplot(d) +
  geom_point(aes(x = sample_loc, y = d13C), size = 3) +
  labs(x = "Catch Location", y = expression(delta~" "^13~C)) +
  theme_classic()

ggplot(d) +
  geom_point(aes(x = sample_loc, y = d15N), size = 3) +
  labs(x = "Catch Location", y = expression(delta~" "^15~N)) +
  theme_classic()

####Gut Contents NEW####
#create a plot showing proprotion of prey items from each benthic habitat type
gut <- read.csv("Lit_Stom_Guts.csv")

#do a lil' data cleaning
gut$Habitat2[gut$Habitat2 == "Rocky/Sandy"] <- "Sandy/Rocky"
gut$Habitat2[gut$Habitat2 == "Rocky/Sandy/OpenWater"] <- "Sandy/Rocky"

gut <- gut[!is.na(gut$Habitat2),]

#summarize the data based on habitat type of prey
gut2 <- gut %>% 
  group_by(Habitat2) %>% 
  summarize(Num = n()) %>% 
  mutate(GutContents = 1)

gut2$Habitat2 <- factor(gut2$Habitat2, levels = c("Open Water", "Rocky", 
                                                "Sandy/Rocky", "Sandy"))

colors <- c("deepskyblue4", "gray27", "khaki4", "peachpuff3")


prey <- ggplot(gut2) +
  geom_bar(aes(x = GutContents, y = Num, fill = Habitat2), 
           position = "fill", stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(y = "Proportion of Prey", x = " ", fill = "Prey Habitat Type") +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white"),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14))

png(filename="MSFigures/PreyHabitat.png", 
    units="in", 
    width=4, 
    height=6, 
    pointsize=8, 
    res=400)

prey

dev.off()

#gray version
colors <- c("deepskyblue4", "gray27", "khaki4", "peachpuff3")


prey_gry <- ggplot(gut2) +
  geom_bar(aes(x = GutContents, y = Num, fill = Habitat2), 
           position = "fill", stat = "identity") +
  scale_fill_grey() +
  labs(y = "Proportion of Prey", x = " ", fill = "Prey Habitat Type") +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white"),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14))

png(filename="MSFigures/PreyHabitat_gry.png", 
    units="in", 
    width=4, 
    height=6, 
    pointsize=8, 
    res=400)

prey_gry

dev.off()

####diversity index?####
library(vegan)
library(BiodiversityR)

d2 <- read.csv("DiversityData.csv")

#removing PDF because that's not actually a species
d2 <- select(d2, !PDF)

# create the species accumulation result
accum <- specaccum(d2, permutations = 100, xlab = "No. of Stomach Samples")

#convert it to a plottable dataframe
accum.long <- as.data.frame(cbind(accum$sites, accum$richness, accum$sd))
names(accum.long) <- c("Stomachs", "Richness", "SD")

#plot it
psa <- ggplot(accum.long, aes(x = Stomachs, y = Richness)) +
  geom_point() +
  stat_smooth(color = "black") + 
  geom_errorbar(aes(ymin=Richness-SD, ymax=Richness + SD), 
                colour="black", width=.4) +
  labs(x = "No. of Stomach Samples", y = "Cumulative No. of Prey Groups") +
  theme_classic()

png(filename="MSFigures/SpecAccum.png", 
    units="in", 
    width=7, 
    height=6, 
    pointsize=8, 
    res=400)

psa

dev.off()


####Gut contents OLD####
gut <- read.csv("Guts_GSBLength.csv", header = T)

#convert sclass to a factor with new labels and arrange species by TL
gut$sclass <- factor(gut$sclass, levels = c("Small", "Medium", "Large"),
                     labels = c("Immature", "Transition", "Mature"))

gut$GutSp <- factor(gut$GutSp, 
                    levels = unique(gut$GutSp[order(gut$preyTL, decreasing = T)]),
                    ordered = T)


#plot by size class with species arranged by TL
ggplot(gut, aes(x = sclass, y = GutSp)) +
  geom_point()

#create dataframe with label values I want later
dlab <- gut %>% 
  select(ID, sclass) %>% 
  distinct() %>% 
  group_by(sclass) %>%
  summarize(value = n())

gut2 <- gut %>% 
  group_by(sclass, GutSp, preyTL) %>% 
  summarize(value = sum(GutNumber))

gut3 <- gut %>% 
  group_by(sclass, GutGroup) %>% 
  summarize(value = sum(GutNumber, na.rm = T),
            meanTL = mean(preyTL))

gut4 <- gut %>%
  group_by(BottomType) %>% 
  summarize(n = sum(GutNumber, na.rm = T))

#convert sclass to a factor with new labels
dlab$sclass <- factor(dlab$sclass, levels = c("Small", "Medium", "Large"),
                      labels = c("Immature", "Transition", "Mature"))
gut2$sclass <- factor(gut2$sclass, levels = c("Small", "Medium", "Large"),
                      labels = c("Immature", "Transition", "Mature"))

gut2$GutSp <- factor(gut2$GutSp, 
                     levels = unique(gut2$GutSp[order(gut2$preyTL, decreasing = T)]),
                     ordered = T)

#plot it!
#create bar label info
ypos <- c(sum(gut2[gut2$sclass == "Mature",]$value, na.rm = T),
          sum(gut2[gut2$sclass == "Transition",]$value, na.rm = T),
          sum(gut2[gut2$sclass == "Immature",]$value, na.rm = T)
)


ggplot(gut2) +
  geom_col(aes(fill = GutSp, y = value, x = sclass)) +
  geom_text(aes(x = sclass, y = ypos, label = value), 
            data = dlab, nudge_y = 1,
            size = 3) +
  labs(x = "Age Group", y = "Count", fill = "Gut Specimen Type\n(Ordered By Trophic Level)") +
  scale_fill_viridis_d() +
  theme_classic()  

#convert sclass to a factor with new labels
gut3$sclass <- factor(gut3$sclass, levels = c("Small", "Medium", "Large"),
                      labels = c("Immature", "Transition", "Mature"))

ypos2 <- c(sum(gut3[gut3$sclass == "Mature",]$value, na.rm = T),
           sum(gut3[gut3$sclass == "Transition",]$value, na.rm = T),
           sum(gut3[gut3$sclass == "Immature",]$value, na.rm = T)
)

gut3$GutGroup <- factor(gut3$GutGroup, 
                        levels = unique(gut3$GutGroup[order(gut3$meanTL, decreasing = T)]),
                        ordered = T)


ggplot(gut3) + 
  geom_col(aes(fill = GutGroup, y = value, x = sclass), 
           position = "stack") +
  geom_text(aes(x = sclass, y = ypos2, label = value), 
            data = dlab, nudge_y = 1,
            size = 3) +
  labs(x = "Age Group", y = "Count", fill = "Gut Specimen Group\n(Ordered by Trophic Level)") +
  scale_fill_viridis_d() +
  theme_classic()

ggplot(gut4) +
  geom_col(aes(y = n, x = BottomType)) +
  labs(x = "Bottom Type", y = "# of Gut Specimens") +
  theme_classic()


ggplot(gut, aes(x = TL_MM, y = preyTL)) +
  geom_point() +
  stat_smooth(method = "lm")
#not looking great.. probs just going to include a table?

ggplot(gut, aes(x = sclass, y = preyTL)) +
  geom_point()