# process the file from the workshop to get the subset of data we want to use
# for the analysis in the paper

# this file is not to be included in the supp. materials

library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(patchwork)


library(sf)
library(stars)
library(units)
library(terra)

load("kit_SST_sandeel.RData")


# subset the data for Shetland
subs <- c("Sumburgh Head", "Noss", "Foula", "Hermaness", "Westerwick", "SE Yell (inc. Burravoe)", "Compass Head", "Noness", "Whale Wick to Sandwick (incl. Ramna Geo)")

shet <- subset(kit, Site %in% subs)
shet$Site <- sub("(\\(.*\\))", "", as.character(shet$Site))
shet$Site <- as.factor(shet$Site)



# out of sample prediction -- Orkney

orkn <- c("Gultak", "Marwick Head", "Mull Head - cliff-nesters",
         "North Hill RSPB, Papa Westray", "Row Head", "Costa Head",
         "Noup Cliffs RSPB (West Westray 2)", "Fair Isle")

ork <- subset(kit, Site %in% orkn)
ork$Site <- sub("(\\(.*\\))", "", as.character(ork$Site))
ork$Site <- as.factor(ork$Site)


save(shet, ork, file="kittiwake_shetland_orkney.RData")
