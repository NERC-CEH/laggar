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


# function to create doughnuts
doughnut_builder <- function(poly_list) {
  # Validate input
  stopifnot(is.list(poly_list), all(map_lgl(poly_list, ~ inherits(.x, "sf"))))
  n <- length(poly_list)
  if (n < 2) stop("Need at least two sf data frames.")
  # First layer
  # results <- list(poly_list[[1]])

  tic("apply")
  diffs <- lapply(1:length(poly_list), function(x) {

    if(i==1){
      return(poly_list)
    } else {
      prev_geoms <- do.call(rbind, poly_list[1:(i - 1)])
      return(st_difference(poly_list[[i]], prev_geoms))
    }

  })
  toc()

  tic("loop")
  diff_outs <- list(poly_list[[1]])

  for(i in 2:length(poly_list)) {
    print(i)
    prev_geoms <- do.call(c, map(poly_list[1:(i - 1)], st_geometry))
    diff_outs[[i]] <- st_difference(poly_list[[i]], prev_geoms)

  }
  toc()

  tic()
  lapply(1:length(poly_list), function(x) {
    if(i==1){
      return(poly_list)
    } else {
      prev_geoms <- do.call(rbind, poly_list[1:(x - 1)])
      geom_diff <- st_difference(poly_list[[x]], prev_geoms)
      geom_diff
    }})
toc()


  the_doughnuts <- map2(poly_list[-1], 2:length(poly_list), function(current_sf, k) {
    # Combine all previous layers' geometries into a single sfc
    prev_geoms <- do.call(c, map(poly_list[1:(k - 1)], st_geometry))
    prev_union <- st_union(st_sfc(prev_geoms, crs = st_crs(current_sf)))

    # Subtract union of previous layers from each polygon in current layer
    geom_diff <- st_difference(st_geometry(current_sf), prev_union)

    # Return sf with same attributes, plus layer ID
    current_sf %>%
      mutate(geometry = geom_diff)

    # # Element-wise st_difference to ensure one geometry per row
    # geom_diff <- purrr::map(
    #   sf::st_geometry(current_sf),
    #   ~ sf::st_difference(sf::st_set_crs(sf::st_sfc(.x), sf::st_crs(current_sf)),
    #                       prev_union)
    # ) |> sf::st_sfc(crs = sf::st_crs(current_sf))
    #
    # # Return sf with updated geometry
    # current_sf$geometry <- geom_diff
    # current_sf
  })
  results <- c(results, the_doughnuts)

  # Combine all results
  results <- do.call(rbind, results)
  results$ID <- 1:nrow(results)
  results
}



# results <- doughnut_builder(poly_list)
#
# ## write code to check doughnut builder
# length(poly_list)
# nrow(results)
#
#
# # check a single object
# nentries <-  unique(sapply(poly_list, nrow))
#
# c(1, 101, 201)
# 1+100
#
#
# rep(1, length(poly_list))
#
#
# ### testing
# doughnut_out = doughnuts
# nbuffers = length(dist)
# object_n = 30

# function to plot single object across all buffer sizes
check_doughnuts <- function(doughnut_out, nbuffers, object_n = 1) {
  if(nbuffers>9)
    stop("Maximum number of buffers for plotting is 9.")

  # get the index
  index_for_plotting <- (0:(nbuffers-1)) *
    (dim(doughnut_out)[1]/nbuffers) + object_n

  # create colour palette
  col_palette <- palette.colors(n = nbuffers)


  # base plot plotting solution
  par(mfrow = c(1, nbuffers+1),
      mai = c(1, 0, 1, 0))
  plot(st_geometry(doughnut_out),
       col = rep(col_palette, each = dim(doughnut_out)[1]/nbuffers),
       main = "All buffered regions")
  lapply(1:length(index_for_plotting), function(x)
    (plot(st_geometry(doughnut_out[index_for_plotting[x],]),
          xlim = st_bbox(doughnut_out[max(index_for_plotting),])[c(1, 3)],
          ylim = st_bbox(doughnut_out[max(index_for_plotting),])[c(2, 4)],
          col = col_palette[x],
          main = paste("Buffered region", x))))
  par(mfrow = c(1, 1))


}


# Function to extract data within an area
# this works with a custom function
## testing
x = seedtraps
y = renv
# if(bng){
dist = c(10, 20, 30, 60, 80, 90, 200)#, 60, 80, 100)

poly_list = buff





doughnut_builder <- function(poly_list) {
  stopifnot(is.list(poly_list), all(purrr::map_lgl(poly_list, ~ inherits(.x, "sf"))))
  n <- length(poly_list)
  if (n < 2) stop("Need at least two sf data frames.")

  # Ensure all layers use the same CRS
  crs_ref <- sf::st_crs(poly_list[[1]])
  poly_list <- purrr::map(poly_list, ~ sf::st_transform(.x, crs_ref))

  results <- list(poly_list[[1]])

  the_doughnuts <- purrr::map2(poly_list[-1], 2:n, function(current_sf, k) {
    prev_geoms <- do.call(c, purrr::map(poly_list[1:(k - 1)], sf::st_geometry))
    prev_union <- sf::st_union(prev_geoms)
    sf::st_crs(prev_union) <- sf::st_crs(current_sf)  # 🔹 enforce CRS

    geom_diff <- purrr::map(
      sf::st_geometry(current_sf),
      ~ sf::st_difference(.x, prev_union)
    ) |> sf::st_sfc(crs = crs_ref)

    current_sf$geometry <- geom_diff
    current_sf
  })

  results <- c(results, the_doughnuts)
  results <- do.call(rbind, results)
  results$ID <- seq_len(nrow(results))
  results
}

#} else {
#   dist = c(10000, 20000, 30000)
# }
func = "sum"


xbuff <- st_buffer(seedtraps[91,], 10)

ggplot() +
  geom_spatraster(data = renv) +
  geom_point(data = seedtraps, aes(x, y), colour = "black") +
  geom_sf(data = xbuff) +
  scale_fill_viridis_c()


# function to do a calculation based on areas within doughnuts
values_within_dist_rast <- function(x, y, dist, func = NULL,
                                    plot_doughnuts = TRUE, object_n = 30) {

  tic("buffer")
  buff <- lapply(dist, sf::st_buffer, x = x)
  toc()

  tic("create doughnuts")
  doughnuts <- doughnut_builder(buff)
  toc()

  if(plot_doughnuts){
    check_doughnuts(doughnuts, nbuffers = length(dist), object_n = object_n)

  }

  if(is.null(func)) {
    tic("extract raw")
    # get values
    #this handles edges i.e.won't count cells that aren't there
    val <- exactextractr::exact_extract(x = y,
                                        y = doughnuts,
                                        fun = func,
                                        include_xy = TRUE,
                                        include_area = TRUE,
                                        include_cols = c("ID"),
                                        include_cell = FALSE)
    toc()

    val <- do.call(rbind, val)

    tic("calculate area")
    # calculate the area of the cell that is covered by the ppolygons
    val$amount_cell_in <- val$cellArea*val$coverage_fraction

    toc()

    return(list(value = val$value, area = val$area*val$coverage_fraction))

  } else if(inherits(func, "character") | inherits(func, "function")){

    tic("extract summary")
    # get values
    #this handles edges i.e.won't count cells that aren't there
    val <- exactextractr::exact_extract(x = y,
                                        y = doughnuts,
                                        fun = func,
                                        force_df = FALSE)
    #convert to matrix
    val <- matrix(val, nrow = nrow(x))
    toc()

    tic("extract area of polys")
    # Calculate area of the doughnuts
    area <- terra::cellSize(y, unit = "m")
    #this handles edges i.e.won't count cells that aren't there
    areaval <- exactextractr::exact_extract(x = area,
                                            y = doughnuts,
                                            fun = "sum",
                                            force_df = FALSE)
    areaval <- matrix(areaval, nrow = nrow(x))
    toc()

    return(list(value = val, area = areaval))

  } else stop("'func' must be of class 'character' or 'function' which takes a vector and outputs a single summary statistic.
         For a full list of available functions see '?exactextractr::exact_extract'.")

}


.sf_conversion <- function(object, name = "object"){
  if(!(any(c("sf","sfc","sfg") %in% class(object)))){
    message(paste("Converting",name,"to sf object"))
    object <- sf::st_as_sf(object)
  }
  return(object)
}

# x = seedtraps
# y = renv
# dist_seq = c(20, 40)
# mindist = NULL
# maxdist = NULL
# incdist = NULL
# bound_poly = NULL

# function to
lg_rast <- function(x, y, dist_seq = NULL, func = NULL, bound_poly = NULL, mindist = NULL,
                    maxdist = NULL, incdist = NULL) {

  # convert x to sf
  x <- .sf_conversion(x, "x")

  # what tests are needed for y ?!?!?!
  if(!inherits(y, "SpatRaster"))
    stop("y must be of class 'SpatRaster'")

  # if (!all(sf::st_geometry_type(x) == "POINT"))
  #   stop("x must be a sf object including only POINT geometries")
  # if (!all(sf::st_geometry_type(y) == "POINT"))
  #   stop("y must be a sf object including only POINT geometries")
  # if (!is.null(bound_poly)) {
  #   bound_poly <- .sf_conversion(bound_poly, "bound_poly")
  #   if (length(bound_poly) != 1 | length(sf::st_geometry_type(bound_poly)) !=
  #       1)
  #     stop("bound_poly must be a sf object comprising a single polygon")
  #   if (sf::st_geometry_type(bound_poly) != "POLYGON")
  #     stop("bound_poly must be a sf object comprising a single polygon")
  # }

  dist <- .index_check(minindex = mindist, maxindex = maxdist,
                       incindex = incdist, index_seq = dist_seq)

  # calculate values within distance
  sumdist <- values_within_dist_rast(x=x, y=y, func=func, dist = dist)

  # calculate value per area
  value_parea <-sumdist$value/sumdist$area


  return(list(buffer_values = sumdist$value, buffer_area = sumdist$area,
              value_parea = value_parea))

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


