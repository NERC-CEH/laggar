
#' Build doughnut geometries from a list of sf data frames
#'
#' This function takes a list of `sf` data frames representing progressively buffered geometries
#' and returns a new `sf` object where each row contains a doughnut-shaped geometry
#' (i.e., the difference between successive buffers). This
#'
#' @param poly_list A list of `sf` data frames. Each data frame must have the same number of rows,
#'   and represent a different buffer level of the same set of geometries.
#'
#' @return An `sf` data frame containing doughnut geometries for each input feature.
#'   The first geometry in each row is preserved, and the rest are doughnut-shaped differences.
#'
#' @export
doughnut_builder <- function(poly_list) {

  # Validate input
  stopifnot(
    is.list(poly_list),
    all(purrr::map_lgl(poly_list, ~ inherits(.x, "sf"))))

  n <- length(poly_list)
  if (n < 2) stop("Need at least two sf data frames.")

  nrows <- unique(sapply(poly_list, nrow))
  if(length(nrows)>1) stop("Number of rows within each item in poly_list must be the same")


  # Combine all sf objects into one
  combined <- do.call(rbind, poly_list)

  # Indexing for each feature across buffer levels
  object_indices <- lapply(1:nrows,
                           function(x) x + (0:(length(poly_list)-1) * nrows))

  # calculate the differences - create doughnuts
  the_doughnuts <- do.call(rbind, purrr::map(object_indices, function(idx) {
    single_object <- combined[idx, ]

    # extract geometries only
    geoms <- sf::st_geometry(single_object)

    # Compute doughnuts: outer buffer minus inner buffer
    diffs <- purrr::map2(
      geoms[-1],  # outer buffers
      geoms[-length(geoms)],  # inner buffers
      sf::st_difference,
      by_feature = TRUE
    )

    # give the differences the same crs and combine
    geom_col <- do.call(c, lapply(diffs, sf::st_sfc, crs = sf::st_crs(geoms)))

    # add the initial object
    sf::st_geometry(single_object) <- c(geoms[1,], geom_col)
    single_object

  }))

  # Combine results and assign IDs
  results <- the_doughnuts[order(as.numeric(row.names(the_doughnuts))), ]
  results$ID <- 1:nrow(results)
  return(results)
}



#' Plot all the buffered polygons and a single object
#'
#' @param doughnut_out Takes the output from `doughnut_builder`. An `sf` data frame.
#' @param nbuffers The number of buffers implemented.
#' @param dist_seq Numeric vector. Distances (in map units) for buffering each measurement point.
#'   Each value represents the outer radius of a doughnut. Must be unique and positive.
#' @param nobjects The number of measurement points that were originally buffered.
#' @param object_n The measurement point that you want to plot separately. A single or vector of integers.
#' @param bg_layer Background layer to plot along with the objects. Currently only accepts `terra::spatRaster`.
#' @return Plots showing all buffered objects and the different buffers for the object chosen in `object_n`.
#' Returns list containing ggplot objects
#' @export
doughnut_checker <- function(doughnut_out, # output of doughnut_builder
                             nbuffers, # number of buffers implemented
                             dist_seq, # buffer distances
                             nobjects, # number of objects/measurement points that were originally buffered
                             object_n = 1, # the number of object that you want to check (i.e. which buffered point do you want plotted). Can be NULL if want to skip single object check, not recommended. Currently might be coded for plotting > 1 object.
                             bg_layer = NULL) # background raster layer for more plotting context
{

  if(nbuffers>9)
    message("Maximum number of buffers for plotting is 9. Recycling colours, plot may look odd.")
  if(max(object_n) > nobjects | any(object_n > nobjects ))
    stop("The object chosen for plotting must be in range 0:nobjects")

  # create an index for plotting
  index_for_plotting <- lapply(object_n, function(x) (0:(nbuffers-1) * nobjects) + x)
  message("plotting object indices ", paste(index_for_plotting, collapse = ", "))

  # create colour palette and add to plotting dataframe
  col_palette <- palette.colors(n = nbuffers, recycle = TRUE)
  colpal <- rep(col_palette, each = dim(doughnut_out)[1]/nbuffers)
  doughnut_out$buffer_name <- rep(factor(dist_seq), each = dim(doughnut_out)[1]/nbuffers)
  names(colpal) <- doughnut_out$buffer_name

  # process raster if supplied
  if(!is.null(bg_layer)){
    if(!inherits(bg_layer, "SpatRaster"))
      stop("Currently not coded for bg_layer to be anything other than 'SpatRaster'")
    if(terra::nlyr(bg_layer) > 1)
      warning("'bg_layer' has more than one layer, using first layer for plotting")
    bg_layer <- as.data.frame(bg_layer,xy = TRUE)
    names(bg_layer) <- c("x", "y", "val")
  }

  # plot all buffers in one plot
  allplot <- ggplot() +
    {if(!is.null(bg_layer)) geom_tile(data = bg_layer, aes(x, y, fill = val))} +
    scale_fill_viridis_c(option = "D", name = NULL) +
    ggnewscale::new_scale_fill() +
    geom_sf(data = doughnut_out, aes(fill =  buffer_name), alpha = 0.7) +
    scale_fill_manual(values =  colpal, name = "Buffer distance") +
    theme_bw()

  ## plot specific layer
  plts <- lapply(1:length(index_for_plotting), function(x)
    lapply(1:length(index_for_plotting[[x]]), function(i)
      ggplot() +
        {if(!is.null(bg_layer)) geom_tile(data = bg_layer, aes(x, y, fill = val))} +
        scale_fill_viridis_c(option = "D", name = NULL) +
        ggnewscale::new_scale_fill() +
        geom_sf(data = doughnut_out[index_for_plotting[[x]][i],],
                aes(fill =  buffer_name), alpha = 0.7) +

        scale_fill_manual(values = colpal, name = NULL, guide = "none") +
        labs(title = paste("Object number", index_for_plotting[[x]][i]),
             subtitle = paste("Buffer distance", unique(doughnut_out[index_for_plotting[[x]][i],]$buffer_name))) +
        coord_sf(xlim = sf::st_bbox(doughnut_out[max(index_for_plotting[[x]]),])[c(1, 3)],
                 ylim = sf::st_bbox(doughnut_out[max(index_for_plotting[[x]]),])[c(2, 4)]) +
        theme_bw()
    )
  )
  obj_plot <- patchwork::wrap_plots(unlist(plts), guides = "collect")
  print(allplot)
  print(obj_plot)

  return(list(allbuffers = allplot,
              objectplots = obj_plot))

}



#' Extract raster values or summaries within concentric buffers ("doughnuts")
#'
#' This function creates concentric buffer zones (doughnuts) around measurement points
#' and extracts raster values within each zone. It can return either raw cell values
#' or summary statistics (e.g., mean, sum) for each doughnut. Optionally, it can plot
#' the doughnuts for visual verification. This is highly recommended.
#'
#' @param x `sf` object. Measurement points around which buffers will be created.
#' @param y Raster object (`terra` SpatRaster). The raster to summarize.
#' @param dist_seq Numeric vector. Distances (in map units) for buffering each measurement point.
#'   Each value represents the outer radius of a doughnut. Must be unique and positive.
#' @param func Summary function applied to raster values within each doughnut.
#'   Default is `NULL`, which returns raw cell values. Can be:
#'   \itemize{
#'     \item A character string naming a function supported by \code{exactextractr::exact_extract}
#'           (e.g., `"mean"`, `"sum"`, `"median"`).
#'     \item A custom function that takes a numeric vector and returns a single value.
#'   }
#' @param plot_doughnuts Logical. If `TRUE`, plots the doughnuts for visual verification.
#'   Recommended to ensure buffering behaves as expected.
#' @param object_n Integer. Row number of the measurement point to plot when
#'   \code{plot_doughnuts = TRUE}. Only a single integer is accepted.
#' @param bg_layer Background layer to plot along with the objects. Currently only accepts `terra::spatRaster`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{value}}{Matrix or vector of extracted values or summary statistics.
#'       Rows correspond to measurement points; columns correspond to doughnuts.}
#'     \item{\code{area}}{Matrix or vector of areas (in square meters) for each doughnut,
#'       accounting for partial cell coverage.}
#'   }
#'
#' @details
#' The function uses \code{exactextractr::exact_extract} to handle partial cell coverage
#' at polygon edges. When \code{func = NULL}, raw cell values are returned along with
#' their corresponding areas. When \code{func} is provided, summary statistics are computed
#' for each doughnut.
#'
#' @examples
#' \dontrun{
#' # Example: Summarize mean raster values within 100m and 500m doughnuts
#' result <- values_within_dist_rast(
#'   x = points_sf,
#'   y = raster_layer,
#'   dist_seq = c(100, 500),
#'   func = "mean",
#'   plot_doughnuts = TRUE
#' )
#'
#' # Access results
#' result$value   # matrix of mean values per doughnut
#' result$area    # matrix of areas per doughnut
#' }
#'
#' @seealso \code{\link[exactextractr]{exact_extract}}, \code{\link[sf]{st_buffer}}, \code{\link[terra]{cellSize}}
#' @export
values_within_dist_rast <- function(x, y, dist_seq, func = NULL,
                                    plot_doughnuts = TRUE, object_n = 1,
                                    bg_layer = NULL) {

  if(any(duplicated(dist_seq))) stop("Cannot have duplicates in 'dist_seq' parameter")

  # buffer sf object
  buff <- lapply(dist_seq, sf::st_buffer, x = x)

  #build the doughnuts
  doughnuts <- doughnut_builder(poly_list = buff)

  # plot the doughnuts to check
  if(plot_doughnuts){
    nobjects <- unique(sapply(buff, nrow))
    if(length(nobjects) > 1)
      stop("cannot have different numbers of objects in each buffer set")
    doughnut_checker(doughnut_out = doughnuts, nbuffers = length(dist_seq), dist_seq = dist_seq,
                     object_n = object_n, nobjects = nobjects, bg_layer = bg_layer)
  }

  if(is.null(func)) {

    # extract raw values
    # this handles edges i.e.won't count cells that aren't there
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

    # extract summary values by polygon
    # this handles edges i.e.won't count cells that aren't there
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



#' Summarize raster within incremental buffer zones around points
#'
#' This function calculates raster values within concentric buffer zones ("doughnuts")
#' around a set of measurement points. It optionally applies a summary function to
#' raster values within each buffer and returns raw values or aggregated summaries.
#'
#' @param x An object representing measurement points. Can be an `sf` object or
#'   convertible to `sf` via internal helper functions.
#' @param y A `SpatRaster` object from the `terra` package representing the raster
#'   to summarize.
#' @param func Optional. A function or character string specifying how to summarize
#'   raster values within each buffer. If `NULL`, raw cell values are returned.
#'   Accepts functions supported by `exactextractr::exact_extract` or a custom
#'   function that takes a numeric vector and returns a single value.
#' @param dist_seq Optional numeric vector of buffer distances. If provided, overrides
#'   `mindist`, `maxdist`, and `incdist`.
#' @param mindist, maxdist, incdist Numeric values specifying the minimum distance,
#'   maximum distance, and increment for generating buffer distances. Ignored if
#'   `dist_seq` is supplied.
#' @param plot_doughnuts Logical. If `TRUE`, plots the doughnut geometries for
#'   visual verification. Highly recommended for error checking.
#' @param bg_layer Background layer to plot along with the objects. Currently only accepts `terra::spatRaster`.
#'
#' @details Internally, this function:
#'   \enumerate{
#'     \item Converts `x` to an `sf` object.
#'     \item Validates that `y` is a `SpatRaster`.
#'     \item Generates buffer distances using `.index_check()`.
#'     \item Builds doughnut geometries and extracts raster values using
#'           `values_within_dist_rast()`.
#'     \item Computes value per unit area for each buffer.
#'   }
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{buffer_values}{Matrix or vector of raster values (raw or summarized) for each buffer.}
#'     \item{buffer_area}{Matrix or vector of buffer areas corresponding to each value.}
#'     \item{value_parea}{Numeric vector of values normalized by buffer area.}
#'   }
#'
#' @seealso \code{\link{values_within_dist_rast}}, \code{\link{doughnut_builder}}, \code{\link{doughnut_checker}}
#' @examples
#' \dontrun{
#' # Example: Summarize mean raster values within 3 buffers around points
#' lg_rast(x = points_sf, y = raster_layer, func = "mean",
#'         mindist = 100, maxdist = 300, incdist = 100)
#' }
#' @export

lg_rast <- function(x,
                    y,
                    func = NULL,
                    dist_seq = NULL, mindist = NULL, maxdist = NULL, incdist = NULL,
                    plot_doughnuts = TRUE, object_n = 1,
                    bg_layer = NULL) {

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
  sumdist <- values_within_dist_rast(x=x, y=y, func=func, dist_seq = dist,
                                     plot_doughnuts = plot_doughnuts,
                                     bg_layer = bg_layer)

  # calculate value per area
  value_parea <-sumdist$value/sumdist$area


  return(list(buffer_values = sumdist$value, buffer_area = sumdist$area,
              value_parea = value_parea))

}




# ### Helper functions ########################
#
# .sf_conversion <- function(object, name = "object"){
#   if(!(any(c("sf","sfc","sfg") %in% class(object)))){
#     message(paste("Converting",name,"to sf object"))
#     object <- sf::st_as_sf(object)
#   }
#   return(object)
# }
#
#
# .index_check <- function(minindex,maxindex,incindex,index_seq){
#   if(any(c(is.null(minindex),is.null(maxindex),is.null(incindex))) & is.null(index_seq))
#     stop("Either index_seq or minindex, maxindex and incindex must be specified")
#   if(is.null(index_seq)){
#     if(!all(is.numeric(minindex),is.numeric(maxindex),is.numeric(incindex)))
#       stop("minindex, maxindex and incindex must be numeric")
#     if(!all(length(minindex)==1,length(maxindex)==1,length(incindex)==1))
#       stop("minindex, maxindex and incindex must be of length 1")
#     if(minindex <= 0){
#       stop("minindex must be greater than 0")
#     }
#     if(maxindex < minindex){
#       stop("maxindex must be greater than minindex")
#     }
#     index <- seq(minindex,maxindex,incindex)
#   } else{
#     if(!is.numeric(index_seq))
#       stop("index_seq must be numeric")
#     if(length(index_seq) < 2)
#       stop("index_seq must be of length greater than 1")
#     if(min(index_seq)<= 0){
#       stop("Smallest indexance increment must be greater than 0")
#     }
#     index <- index_seq
#   }
#   return(index)
# }
