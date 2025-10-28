rm(list = ls())

library(terra)
library(prioritizr)
library(tidyterra)
library(laggar)
library(sf)
library(ggplot2)
library(mgcv)
library(purrr)
library(tictoc)

## choose projection, BNG or latlong
bng = FALSE

## simulate data --- metres
set.seed(100230)
# if(bng){
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


# } else {
#   ## simulate data --- latlong
#   set.seed(100230)
#   boundary <- st_polygon(list(matrix(c(0,100,100,0,0,
#                                        0,0,100,100,0),ncol=2)))
#
#   adults <- st_sample(boundary, size = 300)
#
#   seedtraps <- data.frame(
#     x = rep(seq(5,95,10), 10),
#     y = rep(seq(5,95,10), each = 10)
#   ) |>
#     st_as_sf(coords = c("x","y"), remove = FALSE,
#              crs = "EPSG:32630")
#
#   r <- rast(vect(boundary), nrow = 100, ncol = 100, vals = 1, crs = "EPSG:32630")
#
# }

### simulate data ---------------------------------------
# create empty raster

# simulate rainfall data
renv <- simulate_data(x = r, n = 1, scale = 0.8, intensity = 22, sd = 1)
# plot(renv)
# ---------------------------

## plot
ggplot() +
  geom_spatraster(data = renv) +
  # geom_spatvector(data = terra::vect(buff), alpha = 0.5) +
  geom_point(data = seedtraps, aes(x, y), colour = "black") +
  scale_fill_viridis_c()

# function to create doughnuts
doughnut_builder <- function(poly_list) {
  # Validate input
  stopifnot(is.list(poly_list), all(map_lgl(poly_list, ~ inherits(.x, "sf"))))
  n <- length(poly_list)
  if (n < 2) stop("Need at least two sf data frames.")
  # First layer
  results <- list(poly_list[[1]])

  the_doughnuts <- map2(poly_list[-1], 2:n, function(current_sf, k) {
    # Combine all previous layers' geometries into a single sfc
    prev_geoms <- do.call(c, map(poly_list[1:(k - 1)], st_geometry))
    prev_union <- st_union(st_sfc(prev_geoms, crs = st_crs(current_sf)))

    # Subtract union of previous layers from each polygon in current layer
    geom_diff <- st_difference(st_geometry(current_sf), prev_union)

    # Return sf with same attributes, plus layer ID
    current_sf %>%
      mutate(geometry = geom_diff)
  })
  results <- c(results, the_doughnuts)

  # Combine all results
  results <- do.call(rbind, results)
  results$ID <- 1:nrow(results)
  results
}


# Function to extract data within an area
# this works with a custom function
## testing
x = seedtraps
y = renv
# if(bng){
dist = c(10, 20, 30) #} else {
#   dist = c(10000, 20000, 30000)
# }
func = "sum"

#### I'M SO CONFUSED BY THIS BUFFERING....

xbuff <- st_buffer(seedtraps[91,], 10)

ggplot() +
  geom_spatraster(data = renv) +
  # geom_spatvector(data = terra::vect(buff), alpha = 0.5) +
  geom_point(data = seedtraps, aes(x, y), colour = "black") +
  geom_sf(data = xbuff) +
  scale_fill_viridis_c()



sum_within_dist_rast <- function(x, y, dist, func = sum) {

  # tic("buffer")
  # # buffer point
  # buff <- terra::vect(sf::st_buffer(x, dist))
  # toc()

  tic("buffer")
  buff <- lapply(dist, sf::st_buffer, x = x)
  toc()

  # ggplot() +
  #   geom_tile(data = y, aes(x, y, fill = lyr.1)) +
  #   geom_sf(data = buff[[3]], alpha = 1) +
  #   scale_fill_viridis_c()

  tic("create doughnuts")
  ##### DOES THIS NEED A PLOT TO CHECK THAT IT'S CREATED THE DOUGHNUTS PROPERLY?!?!?! PROBABLY!!
  donuts <- doughnut_builder(buff)
  toc()

  plot(st_geometry(donuts))
  plot(st_geometry(donuts)[1,], add = TRUE, col = "red")
  plot(st_geometry(donuts)[3,], add = TRUE, col = "red")
  plot(st_geometry(donuts)[101,], add = TRUE, col = "blue")
  plot(st_geometry(donuts)[201,], add = TRUE, col = "green")
  plot(st_geometry(donuts)[299,], add = TRUE, col = "green")
  plot(st_geometry(donuts)[300,], add = TRUE, col = "green")

  # # tic("calculate area")
  # # # calculate cell area covered
  # # area <- terra::cellSize(y, unit = "m")
  # # toc()
  #
  # # combine environmental and cell area
  # comb <- c(y, area)
  # names(comb) <- c(names(y), "cellArea")



  if(is.null(func)) {
    tic("extract raw")
    # get values
    #this handles edges i.e.won't count cells that aren't there
    val <- exactextractr::exact_extract(x = y,
                                        y = donuts,
                                        fun = NULL,
                                        include_xy = TRUE,
                                        include_area = TRUE,
                                        include_cols = c("ID"),
                                        include_cell = FALSE)
    toc()

    val <- do.call(rbind, val)
    head(val)

    tic("calculate area")
    # calculate the area of the cell that is covered by the ppolygons
    val$amount_cell_in <- val$cellArea*val$coverage_fraction

    summed_val <- as.numeric(tapply(val[,names(y)], val$ID, func))
    area <- as.numeric(tapply(val$amount_cell_in, val$ID, sum))
    # sumperarea <- summed_val / area
    # mean_val <- as.numeric(tapply(val[,names(y)], val$ID, mean))
    toc()

    list(value = val$value, area = val$area*val$coverage_fraction)

  } else if(inherits(func, "character") | inherits(func, "function")){

    area <- terra::cellSize(y, unit = "m")

    tic("extract summary")
    # get values
    #this handles edges i.e.won't count cells that aren't there
    val <- exactextractr::exact_extract(x = y,
                                        y = donuts,
                                        fun = c("coefficient_of_variation"),
                                        force_df = FALSE)
    head(val)
    toc()

    tic("extract area of polys")
    # get values
    #this handles edges i.e.won't count cells that aren't there
    areaval <- exactextractr::exact_extract(x = area,
                                            y = donuts,
                                            fun = "sum",
                                            force_df = FALSE)
    head(areaval)
    toc()

    list(value = val, area = areaval)

  } else stop("'func' must be of class 'character' or 'function' which takes a vector and outputs a single summary statistic.
         For a full list of available functions see '?exactextractr::exact_extract'.")

  head(val)

  val
  #commenting out to be consistent with sum_within_dist
  #, sumperarea = sumperarea)

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


