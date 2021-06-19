#Kayla Blincow
#4/26/2021

#The purpose of this script is to generate the tables for my isotope chapter.


#clear my workspace
rm(list = ls())

#load packages
library(tidyverse)
library(gt)
library(greekLetters)

#####Gut Contents Table####
gut <- read.csv("PreyTableData.csv", header = T)

names(gut) <- c("Confirmation", "Type", "Common Name", "Scientific Name", 
                "Habitat", "GSB TL (cm)", "Source")

gut$ID <- 1:nrow(gut)

guts <- gut %>% 
  group_by(Confirmation) %>% 
  gt(groupname_col = "Confirmation",
     rowname_col = "ID") %>% 
  tab_style(
    style = list(
      cell_fill(color = "black", alpha = 0.2),
      cell_borders(
        side = c("top", "bottom"),
        color = "black"
      )
    ),
    locations = cells_row_groups()
  ) %>% 
  tab_style(
    style = list(
      cell_text(style = "italic")
    ),
    locations = cells_body(
      columns = "Scientific Name"
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(gt::everything())
  ) %>% 
  tab_options(
    table.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width= px(3)
  ) %>% 
  cols_width(
    starts_with("GSB") ~ px(120),
    "Source" ~ px(200)
  )

gtsave(guts, "MSFigures/GutContentTable.png")


####model table####
m <- read.csv("LooTableData.csv")

m2 <- data.frame(lapply(m, function(x) {
                    gsub("?", greeks("delta"), fixed = T, x)
                }))

loo <- m2 %>% 
  group_by(Isotope) %>% 
  gt() %>% 
  tab_style(
    style = list(
      cell_fill(color = "black", alpha = 0.2),
      cell_borders(
        side = c("top", "bottom"),
        color = "black"
      )
    ),
    locations = cells_row_groups()
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(gt::everything())
  ) %>% 
  tab_options(
    table.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width= px(3)
  ) %>% 
  cols_align(
    align = "center",
    columns = 3:5
  )

gtsave(loo, "MSFigures/LooTable.png")


####Supplement Sample Locations####
d <- read.csv("FinalGSBBulk.csv")

d2 <- d %>% group_by(sample_loc) %>% 
  summarize(n = n()) 

d2 <- arrange(d2, desc(n))
names(d2) <- c("Sample Site", "n")

site <- gt(d2)

gtsave(site, "MSFigures/SuppSampleSite.png")

####Supplement PP Literature Sources####
pp <- read.csv("PPTable.csv")

names(pp) <- c("Type", "Sample Site", paste0("Mean ", greeks("delta"),"13C"), 
               paste0("Mean ", greeks("delta"),"15N"), "Source")

pp$group <- ifelse(pp[,1] == "POM", "Phytoplankton", "Kelp")

pptab <- pp %>% 
  group_by(group) %>% 
  gt() %>% 
  cols_label(Type = "") %>% 
  tab_style(
    style = list(
      cell_fill(color = "black", alpha = 0.2),
      cell_borders(
        side = c("top", "bottom"),
        color = "black"
      )
    ),
    locations = cells_row_groups()
  )  %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(gt::everything())
  ) %>% 
  tab_options(
    table.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width= px(3)
  ) %>% 
  cols_width(
    starts_with("Mean") ~ px(100)
  ) %>% 
  cols_align(
    align = "center",
    columns = 3:4
  )

gtsave(pptab, "MSFigures/SuppLitPPs.png")

