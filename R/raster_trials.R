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


# function to do a calculation based on areas within doughnuts
values_within_dist_rast <- function(x, y, dist, func = NULL,
                                    plot_doughnuts = TRUE, object_n = 1) {

  if(any(duplicated(dist))) stop("Cannot have duplicates in 'dist' parameter")

  # buffer sf object
  buff <- lapply(dist, sf::st_buffer, x = x)

  #build the doughnuts
  doughnuts <- doughnut_builder(poly_list = buff)

  # plot the doughnuts to check
  if(plot_doughnuts){
    nobjects <- unique(sapply(buff, nrow))
    if(length(nobjects) > 1)
      stop("cannot have different numbers of objects in each buffer set")
    doughnut_checker(doughnuts, nbuffers = length(dist), object_n = object_n,
                     nobjects = nobjects)
  }

  if(is.null(func)) {

    # extract raw values
    # get values
    #this handles edges i.e.won't count cells that aren't there
    val <- exactextractr::exact_extract(x = y,
                                        y = doughnuts,
                                        fun = func,
                                        include_xy = TRUE,
                                        include_area = TRUE,
                                        include_cols = c("ID"),
                                        include_cell = FALSE)

    val <- do.call(rbind, val)

    return(list(value = val$value, area = val$area*val$coverage_fraction))

  } else if(inherits(func, "character") | inherits(func, "function")){

    # extract summary values
    # get values
    #this handles edges i.e.won't count cells that aren't there
    val <- exactextractr::exact_extract(x = y,
                                        y = doughnuts,
                                        fun = func,
                                        force_df = FALSE)
    #convert to matrix
    val <- matrix(val, nrow = nrow(x))

    # extract area of each doughnut
    # Calculate area of the doughnuts
    area <- terra::cellSize(y, unit = "m")
    #this handles edges i.e.won't count cells that aren't there
    areaval <- exactextractr::exact_extract(x = area,
                                            y = doughnuts,
                                            fun = "sum",
                                            force_df = FALSE)
    areaval <- matrix(areaval, nrow = nrow(x))

    return(list(value = val, area = areaval))

  } else stop("'func' must be of class 'character' or 'function' which takes a vector
  and outputs a single summary statistic.
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
lg_rast <- function(x,
                    y,
                    func = NULL, # function to summarise the data within each of the buffers (as specified in exactextractr::exact_extract. If NULL it returns raw values for each cell within the buffer
                    #bound_poly = NULL, -- no need for bound poly because the area is limited by the raster
                    dist_seq = NULL, mindist = NULL, maxdist = NULL, incdist = NULL, # buffer specification
                    plot_doughnuts = TRUE # plot the doughnuts for error checking - highly recommended!!
) {

  # convert x to sf
  x <- .sf_conversion(x, "x")

  # what tests are needed for x and y ?!?!?!
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
  sumdist <- values_within_dist_rast(x=x, y=y, func=func, dist = dist,
                                     plot_doughnuts = plot_doughnuts)

  # calculate value per area
  value_parea <-sumdist$value/sumdist$area


  return(list(buffer_values = sumdist$value, buffer_area = sumdist$area,
              value_parea = value_parea))

}








