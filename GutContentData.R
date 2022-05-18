#Kayla Blincow
# 4/15/2020
#UPDATE: 4/28/2022 - taking a look at other summaries of the data for revisions

#Data exploration of Mike Shane's GSB gut content data.

#clear my workspace
rm(list = ls())

library(tidyverse)

#load split gut content data
d <- read.csv("gutcontents_split.csv", header = T)

d$Name <- as.factor(d$Name)
d$Family <- as.factor(d$Family)
d$Order <- as.factor(d$Order)

#calculate percentage of unique prey items that are fish versus inverts
d2 <- d %>% 
  group_by(Type) %>% 
  summarize(ntype = n_distinct(Name))

d2$ntype[1]/sum(d2$ntype)
#60% are fish, 40% are invertebrates

#verify habitat type
d3 <- d %>% 
  group_by(Habitat) %>% 
  summarize(ntype = n_distinct(Name, na.rm = T))

1 - d3$ntype[1]/sum(d3$ntype) #87.5% associated with benthic habitats

d3$ntype/sum(d3$ntype)

#look at prey size versus GSB size
ggplot(d, aes(x = PreyLength, y = GSBLength)) +
  geom_point()
#not a whole lot there.. don't have enough data..

#look at relative percentages of different families/orders
d4 <- d %>% 
  group_by(Family) %>% 
  summarize(ntype = n_distinct(Name))

d5 <- d %>% 
  group_by(Type,Order) %>% 
  summarize(ntype = n_distinct(Name))

filter(d5, Type == "Fish" & Order != "NA") %>% 
  summarize(ntype/sum(ntype))

filter(d5, Type == "Invertebrate" & Order != "NA") %>% 
  summarize(ntype/sum(ntype))



#OLD STUFF
#read in the data
d <- read.csv("HSWRIGutsAdapted.csv")
d <- d[,c(1:6,9:15)]
d <- d[-31,]

d$MORT_DATE <- as.Date(d$MORT_DATE, format = "%m/%d/%Y")
d$sizeclass <- NA

for(i in 1:length(d$sizeclass)){
  if(d$TL_MM[i] < 500){
    d$sizeclass[i] <- "0-500"
  } else if(d$TL_MM[i] > 500 & d$TL_MM[i] < 1000) {
    d$sizeclass[i] <- "500-1000" 
  } else if(d$TL_MM[i] > 1000 & d$TL_MM[i] < 1500) {
    d$sizeclass[i] <- "1000-1500"
  } else {
    d$sizeclass[i] <- "1500-2000"
  }
  
}

d$sizeclass <- as.factor(d$sizeclass)

plot(x = d$HabitatGroup, y = d$TL_MM)
plot(x = d$BottomType, y = d$TL_MM)

#what kind of plots do I want to use to visualize this??
ggplot(d) +
  geom_point(aes(x = GutSp, y = TL_MM, shape = GutGroup, color = HabitatGroup),
             size = 2)

d2 <- d %>% 
  group_by(sizeclass, BottomType) %>% 
  summarize(count = n()) %>% 
  drop_na()


d2$GutContents <- as.factor("")

colors <- c("deepskyblue4", "gray27", "peachpuff3", "khaki4")

ggplot(d2) +
  geom_bar(aes(x = GutContents, y = count, fill = BottomType), 
           position = "fill", stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(y = "Percent of Gut Contents", x = " ", fill = "Prey Habitat Type") +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white"),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14))

ggsave("HabtiatPlot.png", units="in", width=4, height=6, dpi=1000)
