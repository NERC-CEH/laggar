rm(list = ls())

library(terra)
library(prioritizr)
library(tidyterra)
library(laggar)
library(sf)
library(ggplot2)
library(mgcv)
library(purrr)

## simulate data
set.seed(100230)
boundary <- st_polygon(list(matrix(c(0,1000,1000,0,0,
                                     0,0,1000,1000,0),ncol=2)))

adults <- st_sample(boundary, size = 300)

# st_crs(boundary) <- "EPSG:27700"

seedtraps <- data.frame(
  x = rep(seq(50,950,100), 10),
  y = rep(seq(50,950,100), each = 10)
) |>
  st_as_sf(coords = c("x","y"), remove = FALSE,
           crs = "EPSG:27700")


# ggplot() + geom_sf(data = boundary) + geom_sf(data = adults) +
#   geom_sf(data = seedtraps, colour = "red") +
#   theme_void()

# rasterize(boundary)


### simulate data ---------------------------------------
# create empty raster
r <- rast(vect(boundary), nrow = 100, ncol = 100, vals = 1, crs = "EPSG:27700")

# simulate rainfall data
renv <- simulate_data(x = r, n = 1, scale = 0.8, intensity = 22, sd = 1)
plot(renv)
# ---------------------------

## plot
ggplot() +
  geom_spatraster(data = renv) +
  # geom_spatvector(data = terra::vect(buff), alpha = 0.5) +
  geom_point(data = seedtraps, aes(x, y), colour = "black") +
  scale_fill_viridis_c()


##### This is a function that creates donuts
doughnut_builder <- function(poly_list) {
  # Ensure input is a list of sfc or sfg objects
  stopifnot(is.list(poly_list))
  n <- length(poly_list)

  if (n < 2) stop("Need at least two polygons.")

  # Apply st_difference sequentially
  diffs <- map2( poly_list[-1], poly_list[-n], st_difference)

  # Return as an sfc geometry collection
  diffs <- do.call(rbind, diffs)

  #remove second coordinates
  diffs[, grep(".1", colnames(diffs))] <- NULL
  diffs
}





# Function to extract data within an area
# this works with a custom function
## testing
x = seedtraps
y = renv
dist = c(10, 20)
func = sum

sum_within_dist_rast <- function(x, y, dist, func = sum) {

  # buffer point
  buff <- terra::vect(sf::st_buffer(x, dist))

  buff <- lapply(dist, sf::st_buffer, x = x)

  # ggplot() +
  #   geom_tile(data = y, aes(x, y, fill = lyr.1)) +
  #   geom_sf(data = buff[[3]], alpha = 1) +
  #   scale_fill_viridis_c()

  donuts <- doughnut_builder(buff)
  buff_donuts <- rbind(buff[[1]], donuts)

  # plot(st_geometry(buff_donuts[5,]), col = c("black"))
  # plot(st_geometry(buff_donuts[4,]), col = c("yellow"), add = TRUE)
  # plot(st_geometry(buff_donuts[3,]), col = c("green", "yellow", "black"), add = TRUE)
  # plot(st_geometry(buff_donuts[2,]), col = c("blue", "green", "yellow", "black"), add = TRUE)
  # plot(st_geometry(buff_donuts[1,]), col = c("red", "blue", "green", "yellow", "black"), add = TRUE)

  # calculate cell area covered
  area <- terra::cellSize(y, unit = "m")

  # combine environmental and cell area
  comb <- c(y, area)

  # get values
  #this handles edges i.e.won't count cells that aren't there
  system.time(val <- terra::extract(comb, buff_donuts,
                                    xy = TRUE,
                                    touches = TRUE, # all cells that are touched by buffer, not just centroid
                                    exact = TRUE))

  # calculate the area of the cell that is covered by the ppolygons
  val$amount_cell_in <- val$area*val$fraction

  summed_val <- as.numeric(tapply(val[,names(y)], val$ID, func))
  area <- as.numeric(tapply(val$amount_cell_in, val$ID, sum))
  # sumperarea <- summed_val / area
  # mean_val <- as.numeric(tapply(val[,names(y)], val$ID, mean))

  # commenting out to be consistent with sum_within_dist
  list(summed_val = summed_val, area = area)#, sumperarea = sumperarea)

}

library(sf)
library(purrr)



# polys <- list(buff)
#
# # Sequential difference
# donuts <- sequential_difference(buff)
#
# plot(st_geometry(donuts), col = c("lightblue", "lightgreen"), border = "black")


############### Original POINTS code uses lg_points
############### I think changing sum_within_dist function is enough
# lg_points <-
# function (x, y, dist_seq = NULL, bound_poly = NULL, mindist = NULL,
#     maxdist = NULL, incdist = NULL)
# {
#     x <- .sf_conversion(x, "x")
#     y <- .sf_conversion(y, "y")
#     if (!all(sf::st_geometry_type(x) == "POINT"))
#         stop("x must be a sf object including only POINT geometries")
#     if (!all(sf::st_geometry_type(y) == "POINT"))
#         stop("y must be a sf object including only POINT geometries")
#     if (!is.null(bound_poly)) {
#         bound_poly <- .sf_conversion(bound_poly, "bound_poly")
#         if (length(bound_poly) != 1 | length(sf::st_geometry_type(bound_poly)) !=
#             1)
#             stop("bound_poly must be a sf object comprising a single polygon")
#         if (sf::st_geometry_type(bound_poly) != "POLYGON")
#             stop("bound_poly must be a sf object comprising a single polygon")
#     }
#     dist <- .index_check(minindex = mindist, maxindex = maxdist,
#         incindex = incdist, index_seq = dist_seq)
#     sumdist <- lapply(dist, sum_within_dist, x = x, y = y, bound_poly = bound_poly)
#     num_points <- sapply(sumdist, "[[", 1)
#     num_points <- t(apply(num_points, 1, function(x) diff(c(0,
#         x))))
#     if (!is.null(bound_poly)) {
#         buffer_area <- sapply(sumdist, "[[", 2)
#         buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,
#             x))))
#         points_parea <- num_points/buffer_area
#     }
#     else {
#         buffer_area <- NA
#         points_parea <- NA
#     }
#     return(list(num_points = num_points, buffer_area = buffer_area,
#         points_parea = points_parea))
# }


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

lg_rast <- function(x, y, dist_seq = NULL, func = "sum", bound_poly = NULL, mindist = NULL,
                    maxdist = NULL, incdist = NULL) {
  x <- .sf_conversion(x, "x")
  y <- .sf_conversion(y, "y")
  if (!all(sf::st_geometry_type(x) == "POINT"))
    stop("x must be a sf object including only POINT geometries")
  if (!all(sf::st_geometry_type(y) == "POINT"))
    stop("y must be a sf object including only POINT geometries")
  if (!is.null(bound_poly)) {
    bound_poly <- .sf_conversion(bound_poly, "bound_poly")
    if (length(bound_poly) != 1 | length(sf::st_geometry_type(bound_poly)) !=
        1)
      stop("bound_poly must be a sf object comprising a single polygon")
    if (sf::st_geometry_type(bound_poly) != "POLYGON")
      stop("bound_poly must be a sf object comprising a single polygon")
  }

  dist <- .index_check(minindex = mindist, maxindex = maxdist,
                       incindex = incdist, index_seq = dist_seq)
  sumdist <- lapply(dist, sum_within_dist_rast, x = x, y = y, func = func)#, bound_poly = bound_poly)

  sum_within_dist_rast(x=x, y=y, func=func, dist = dist_seq)

  num_points <- sapply(sumdist, "[[", 1)
  num_points <- t(apply(num_points, 1, function(x) diff(c(0,
                                                          x))))
  if (!is.null(bound_poly)) {
    buffer_area <- sapply(sumdist, "[[", 2)
    buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,
                                                              x))))
    points_parea <- num_points/buffer_area
  }
  else {
    buffer_area <- NA
    points_parea <- NA
  }
  return(list(num_points = num_points, buffer_area = buffer_area,
              points_parea = points_parea))

}























dist <- laggar:::.index_check(minindex = NULL, maxindex = NULL,
                              incindex = NULL, index_seq = NULL)


sum_within_dist(x = seedtraps, y = adults, poly = poly, dist = 20)


dist_seq <- c(20)

### this works but takes a while with big buffer areas because of extract function presumably
sumdist <- sum_within_dist_rast(x = seedtraps,
                                y = renv,
                                dist = 100)

## sum the number of points within dist of each seedtrap
sumdist <- lapply(dist_seq, sum_within_dist, x = x, y = y, poly = poly)

## no idea what this does
num_points <- sapply(sumdist, "[[", 1) # extracts the first item in each sublist, i.e. npoints
num_points <- t(apply(num_points, 1, function(x) diff(c(0,
                                                        x))))


if (!is.null(poly)) {
  buffer_area <- sapply(sumdist, "[[", 2)
  buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,
                                                            x))))
  points_parea <- num_points/buffer_area
} else {
  buffer_area <- NA
  points_parea <- NA
}
return(list(num_points = num_points, buffer_area = buffer_area,
            points_parea = points_parea))


tst

points_perha <- 1e4*tst$sumperarea

lgr <- lg_assembly(response = nseeds,
                   covariate = points_perha,
                   index_seq = dist_seq)


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


