
library(terra)
library(prioritizr)
library(tidyterra)
# library(laggar)
library(sf)
library(ggplot2)
library(mgcv)
library(purrr)
library(tictoc)

# source("R/raster.R")

## simulate spatial data
set.seed(100230)
boundary <- st_polygon(list(matrix(c(0,1000,1000,0,0,
                                     0,0,1000,1000,0),ncol=2)))

adults <- st_sample(boundary, size = 300)

seedtraps <- data.frame(
  x = rep(seq(50,950,100), 10),
  y = rep(seq(50,950,100), each = 10)
) |>
  st_as_sf(coords = c("x","y"), remove = FALSE,
           crs = "EPSG:32630")

r <- rast(vect(boundary), nrow = 100, ncol = 100, vals = 1, crs = "EPSG:32630")


### simulate data ---------------------------------------
# create empty raster

# simulate rainfall data
renv <- simulate_data(x = r, n = 1, scale = 0.8, intensity = 22, sd = 1)
# ---------------------------

## plot
ggplot() +
  geom_spatraster(data = renv) +
  # geom_spatvector(data = terra::vect(buff), alpha = 0.5) +
  geom_point(data = seedtraps, aes(x, y), colour = "black") +
  scale_fill_viridis_c()

# # testing
# poly_list = buff



# Function to extract data within an area
# this works with a custom function
## testing
x = seedtraps
y = renv
dist = c(10,20,30)
object_n = 1
plot_doughnuts = TRUE
func = "sum"


xbuff <- st_buffer(seedtraps[91,], 10)

ggplot() +
  geom_spatraster(data = renv) +
  geom_point(data = seedtraps, aes(x, y), colour = "black") +
  geom_sf(data = xbuff) +
  scale_fill_viridis_c()




.sf_conversion <- function(object, name = "object"){
  if(!(any(c("sf","sfc","sfg") %in% class(object)))){
    message(paste("Converting",name,"to sf object"))
    object <- sf::st_as_sf(object)
  }
  return(object)
}

x = seedtraps
y = renv
dist_seq = c(20, 40)
mindist = NULL
maxdist = NULL
incdist = NULL
bound_poly = NULL






