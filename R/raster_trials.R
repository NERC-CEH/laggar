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
                                    plot_doughnuts = TRUE, object_n = 30) {

  if(any(duplicated(dist))) stop("Cannot have duplicates in 'dist' parameter")

  tic("buffer")
  buff <- lapply(dist, sf::st_buffer, x = x)
  toc()

  tic("create doughnuts")
  doughnuts <- doughnut_builder(poly_list = buff)
  toc()

  if(plot_doughnuts){
    ### ADD THE POSSIBILITY TO PLOT IT WITH THE ENVIRONMENTAL BACKGROUND
    ### CODE IN GGPLOT
    doughnut_checker(doughnuts, nbuffers = length(dist), object_n = object_n,
                     nobjects = unique(sapply(buff, nrow)))
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






# function to create doughnuts
doughnut_builder <- function(poly_list) {
  # Validate input
  stopifnot(is.list(poly_list), all(map_lgl(poly_list, ~ inherits(.x, "sf"))))
  n <- length(poly_list)
  if (n < 2) stop("Need at least two sf data frames.")

  tic("get differences")

  nrows <- unique(sapply(poly_list, nrow))
  if(length(nrows)>1) stop("Number of rows within each item in poly_list must be the same")

  # combine the list into a single object to work on
  comb_list <- do.call(rbind, poly_list)

  # identify each buffered object as determined by number of polygons supplied
  # need to get the buffered objects of object 1:n objects in poly_list
  object_indices <- lapply(1:nrows, function(x) x + (0:(length(poly_list)-1) * nrows))

  # tic("my code")
  # # get differences between each polygon
  # the_doughnuts <- do.call(rbind, lapply(object_indices, function(x) {
  #   single_object <- comb_list[x, ]
  #   rbind(single_object[1,], do.call(rbind, lapply(1:nrow(single_object), function(i) {
  #
  #     prev_geoms <- st_union(single_object[1:(i - 1),])
  #     st_difference(single_object[i,], prev_geoms)
  #   }))
  #   )
  # }))
  # toc()

  # # code with help from copilot
  # tic("chat code")

  # Assuming each object has n buffers stacked in order

  the_doughnuts <- do.call(rbind, map(object_indices, function(idx) {
    single_object <- comb_list[idx, ]

    # Extract geometries only
    geoms <- st_geometry(single_object)

    # Compute doughnuts: outer buffer minus inner buffer
    diffs <- map2(
      geoms[-1],  # outer buffers
      geoms[-length(geoms)],  # inner buffers
      st_difference,
      by_feature = TRUE
    )

    geom_col <- do.call(c, lapply(diffs, st_sfc, crs = st_crs(geoms)))


    st_geometry(single_object) <- c(geoms[1,], geom_col)
    single_object

  }))

  toc()


  # # })
  #
  # # remove extra cols
  # for(i in 1:length(the_doughnuts)) {
  #   the_doughnuts[[i]][,grep(".1", colnames(the_doughnuts[[i]]), value = FALSE)] <- NULL
  # }

  results <- the_doughnuts[ order(as.numeric(row.names(the_doughnuts))), ] # do.call(rbind, the_doughnuts)

  # # what is going on with this geometry????
  # plot(st_geometry(results[301,]))
  # plot(st_geometry(results[301:400,]), col = "red")
  # plot(st_geometry(results[201:300,]), col = "blue", add = TRUE)
  # plot(st_geometry(results[101:200,]), col = "yellow", add = TRUE)
  # plot(st_geometry(results[1:100,]), col = "green", add = TRUE)

  results$ID <- 1:nrow(results)
  results
}


# ### testing
# doughnut_out = doughnuts
# nbuffers = length(dist)
# nobjects = unique(sapply(buff, nrow))
# object_n = 30

# function to plot single object across all buffer sizes
doughnut_checker <- function(doughnut_out, # output of doughnut_builder
                             nbuffers, # number of buffers implemented
                             nobjects, # number of objects that were originally buffered
                             object_n = 1) # the number of object that you want to check (i.e. which buffered point do you want plotted)
{
  if(nbuffers>9)
    stop("Maximum number of buffers for plotting is 9.")
  if(object_n > nobjects)
    stop("The object chosen for plotting must be in range 0:nobjects")


  # create an index for plotting
  index_for_plotting <-
    (0:(nbuffers-1) * nobjects) + object_n
  message("plotting object indices ", paste(index_for_plotting, collapse = ", "))

  # # get the index
  # index_for_plotting <- (0:(nbuffers-1)) *
  #   (dim(doughnut_out)[1]/nbuffers) + object_n

  # create colour palette
  col_palette <- palette.colors(n = nbuffers)


  # base plot plotting solution
  plot(st_geometry(doughnut_out),
       col = rep(col_palette, each = dim(doughnut_out)[1]/nbuffers),
       main = "All buffered regions")

  par(mfrow = c(1, nbuffers))
  lapply(1:length(index_for_plotting), function(x)
    (plot(st_geometry(doughnut_out[index_for_plotting[x],]),
          xlim = st_bbox(doughnut_out[max(index_for_plotting),])[c(1, 3)],
          ylim = st_bbox(doughnut_out[max(index_for_plotting),])[c(2, 4)],
          col = col_palette[x],
          main = paste("Buffered region", x))))
  par(mfrow = c(1, 1))


}








