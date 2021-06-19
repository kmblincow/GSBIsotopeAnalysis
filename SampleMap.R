#Kayla Blincow
#4/16/2020

#The purpose of this script is create a map of the sampling locations for GSB
#fin clips

#clear my workspace
rm(list = ls())

#load the libraries
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(dplyr)
library(sf)
library(ggsn)

#load the data
d <- read.csv("FinalGSBBulk.csv", header = T)

#find lat/long data
d2 <- select(d, lat, long, sample_loc) %>% 
  unique()

d2 <- d2[-c(2, 10, 14:16),]


#breakdown lat/long and sum
d3 <- d %>% 
  group_by(sample_loc) %>% 
  summarize(n()) %>% 
  left_join(d2)

names(d3) <- c("sample_loc", "tot", "lat", "long")

#make the larger Caribbean map for context
world <- map_data("mapdata::worldHires")
world <- fortify(world)


map <- ggplot() + 
  geom_polygon(data = world, aes(long, lat, group = group),
                               fill = "gray60", color = "gray30") +
  geom_point(data = d3, aes(x = long, y = lat), size = 4, shape = 17) +
  coord_cartesian(xlim = c(-121, -105), ylim = c(25, 35)) +
  north(symbol = 16, location = "topright", scale = 0.2, 
        x.min = -121, x.max = -105, y.min = 25, y.max = 35) +
  scalebar(dist = 250, dist_unit = "km", x.min = -121, 
           x.max = -105, y.min = 25, y.max = 35, location = "bottomleft",
           transform = TRUE, model = "WGS84", st.bottom = FALSE,
           st.dist = 0.03, st.size = 4) +
  labs(x = "Longitude", y = "Latitude")+
  theme_classic()+
  theme(text = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) 


png(filename="MSFigures/Iso_SampleMap.png", 
units="in", 
width=7, 
height=6, 
pointsize=8, 
res=400)

map

dev.off()
