rm(list = ls())

library(terra)
library(prioritizr)
library(tidyterra)
library(laggar)
library(sf)
library(ggplot2)
library(mgcv)

## simulate data
set.seed(100230)
boundary <- st_polygon(list(matrix(c(0,1000,1000,0,0,
                                     0,0,1000,1000,0),ncol=2)))

adults <- st_sample(boundary, size = 300)

seedtraps <- data.frame(
  x = rep(seq(50,950,100), 10),
  y = rep(seq(50,950,100), each = 10)
) |>
  st_as_sf(coords = c("x","y"), remove = FALSE)


ggplot() + geom_sf(data = boundary) + geom_sf(data = adults) +
  geom_sf(data = seedtraps, colour = "red") +
  theme_void()

# rasterize(boundary)

# create empty raster
r <- rast(vect(boundary), nrow = 100, ncol = 100, vals = 1, crs = "PROJ.4")

# simulate rainfall data
renv <- simulate_data(x = r, n = 1, scale = 0.8, intensity = 22, sd = 1)
plot(renv)

class(renv)

ggplot() +
  geom_spatraster(data = renv) +
  # geom_spatvector(data = terra::vect(buff), alpha = 0.5) +
  geom_point(data = seedtraps, aes(x, y), colour = "black") +
  scale_fill_viridis_c()


# Function to extract sum of data within an area
sum_within_dist_rast <- function(x, y, dist) {


  # ## need to create a coordinate system to try and get area I think...
  # terra::crs(y) <- 27700

  # buffer point
  buff <- terra::vect(sf::st_buffer(x, dist))

  # get values
  val <- terra::extract(y, buff, xy = TRUE,
                        exact = TRUE )

  summed_val <- tapply(sumdist[,names(y)], val$ID, sum)
  # ## fraction doesn't give us area! It gives us proportion of each cell that is covered.
  # area <- tapply(val$fraction, val$ID, sum)
  area <- expanse(buff)
  suppressWarnings(area_in_buff <- sf::st_area(sf::st_intersection(buff, poly)))
  terra::intersect(y, buff)

  ## need to get the intersection between each buffered aread and the raster
  ## and then the area

  # # sum the values ----- CHECK TO SEE ABOUT BOUNDARY CASES
  # out_list <- lapply(1:length(unique(val$ID)), function(x) {
  #   list(summed_value = sum(val$lyr.1[val$ID == unique(val$ID)[x]]),
  #        area =  sum(val$fraction[val$ID == unique(val$ID)[x]]))
  # })
  # outs <- lapply(1:length(out_list[[1]]), function(x) sapply(out_list, "[[", x))
  # names(outs) <- names(out_list[[1]])


  summed_val <- as.numeric(tapply(val[,names(y)], val$ID, sum))
  area <- as.numeric(tapply(val$fraction, val$ID, sum))
  sumperarea <- summed_val / area

  list(summed_val = summed_val, area = area, sumperarea = sumperarea)

  # list(npoints = sum_y_in_buff, area = area_in_buff)

}


### this works but takes a while with big buffer areas because of extract function presumably
tst <- sum_within_dist_rast(x = seedtraps,
                            y = renv,
                            dist = 100)

tst
dist_seq <- seq(10,200,10)

## sum the number of points within dist of each seedtrap
system.time(
  sumdist <- lapply(dist_seq, sum_within_dist_rast, x = seedtraps, y = renv)
)

## no idea what this does
num_points <- sapply(sumdist, "[[", 1) # extracts the first item in each sublist, i.e. npoints
# convert from cumulative sum
num_points <- t(apply(num_points, 1, function(x) diff(c(0,
                                                        x))))

buffer_area <- sapply(sumdist, "[[", 2)
buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,x))))
points_parea <- num_points/buffer_area




### for a sequence
dist_seq <- seq(10,200,10)

x = seedtraps
y = renv
dist = 20#seq(2,10, by = 2)


sumdist <- sum_within_dist_rast(x = seedtraps,
                                y = renv,
                                dist = 20)
glimpse(sumdist)


# buffer point
buff <- terra::vect(sf::st_buffer(x, dist))

# get values
val <- terra::extract(y, buff, xy = TRUE,
                      exact = TRUE )

# sum the values
lapply(1:length(unique(val$ID)), function(x) {
  list(summed_value = sum(val$lyr.1[val$ID == unique(val$ID)[x]]),
       area =  sum(val$fraction[val$ID == unique(val$ID)[x]]))
})

# list(npoints = sum_y_in_buff, area = area_in_buff)


## sum the number of points within dist of each seedtrap
sumdist <- lapply(dist_seq, sum_within_dist_rast, x = seedtraps, y = renv)





x <- rasterize(buff, renv, fun="sum")

# get intersection
r <- terra::relate(x = terra::vect(buff),
                   y = renv,
                   "intersects")
r






y_in_buff <- sf::st_contains(buff, y)
sum_y_in_buff <- sapply(y_in_buff, length)
if(!is.null(poly)){
  suppressWarnings(area_in_buff <- sf::st_area(sf::st_intersection(buff, poly)))
} else{
  area_in_buff <- rep(NA, length(sum_y_in_buff))
}
list(npoints = sum_y_in_buff, area = area_in_buff)



sum_within_dist <- function(x, y, poly, dist){
  buff <- sf::st_buffer(x, dist)
  y_in_buff <- sf::st_contains(buff, y)
  sum_y_in_buff <- sapply(y_in_buff, length)
  if(!is.null(poly)){
    suppressWarnings(area_in_buff <- sf::st_area(sf::st_intersection(buff, poly)))
  } else{
    area_in_buff <- rep(NA, length(sum_y_in_buff))
  }
  list(npoints = sum_y_in_buff, area = area_in_buff)
}


