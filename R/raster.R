

#' Build Doughnut Geometries from a List of sf Data Frames
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
#' @importFrom sf st_geometry st_difference st_crs st_sfc
#' @importFrom purrr map map2 map_lgl
#' @export
doughnut_builder <- function(poly_list) {

  # Validate input
  stopifnot(
    is.list(poly_list),
    all(map_lgl(poly_list, ~ inherits(.x, "sf"))))

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
    geoms <- st_geometry(single_object)

    # Compute doughnuts: outer buffer minus inner buffer
    diffs <- map2(
      geoms[-1],  # outer buffers
      geoms[-length(geoms)],  # inner buffers
      st_difference,
      by_feature = TRUE
    )

    # give the differences the same crs and combine
    geom_col <- do.call(c, lapply(diffs, st_sfc, crs = st_crs(geoms)))

    # add the initial object
    st_geometry(single_object) <- c(geoms[1,], geom_col)
    single_object

  }))

  # Combine results and assign IDs
  results <- the_doughnuts[order(as.numeric(row.names(the_doughnuts))), ]
  results$ID <- 1:nrow(results)
  return(results)
}



# function to plot single object across all buffer sizes
doughnut_checker <- function(doughnut_out, # output of doughnut_builder
                             nbuffers, # number of buffers implemented
                             nobjects, # number of objects/measurement points that were originally buffered
                             object_n = 1, # the number of object that you want to check (i.e. which buffered point do you want plotted). Can be NULL if want to skip single object check, not recommended. Currently might be coded for plotting > 1 object.
                             bg_layer = NULL) # background raster layer for more plotting context
{

  if(nbuffers>9)
    stop("Maximum number of buffers for plotting is 9.")
  if(max(object_n) > nobjects | any(object_n > nobjects ))
    stop("The object chosen for plotting must be in range 0:nobjects")


  # create an index for plotting
  index_for_plotting <- lapply(object_n, function(x) (0:(nbuffers-1) * nobjects) + x)
  # (0:(nbuffers-1) * nobjects) + object_n
  message("plotting object indices ", paste(index_for_plotting, collapse = ", "))

  # create colour palette and add to plotting dataframe
  col_palette <- palette.colors(n = nbuffers)
  colpal <- rep(col_palette, each = dim(doughnut_out)[1]/nbuffers)
  doughnut_out$buffer_name <- rep(factor(dist), each = dim(doughnut_out)[1]/nbuffers)
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
  ggplot() +
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
        coord_sf(xlim = st_bbox(doughnut_out[max(index_for_plotting[[x]]),])[c(1, 3)],
                 ylim = st_bbox(doughnut_out[max(index_for_plotting[[x]]),])[c(2, 4)]) +
        theme_bw()
    )
  )

    patchwork::wrap_plots(unlist(plts), guides = "collect", ncol = nbuffers)

}




# function to do a calculation based on areas within doughnuts

#' @param x The measurement points
#' @param y The raster to summarise
#' @param dist The distances by which to buffer each of the measurement points
#' @param func Function by which to summarise the values within each doughnut. Default
#'      is NULL which returns the raw values. Can be a function listed in exactextractr::exact_extract
#'      or a custom function which takes a vector and returns a single value
#' @param plot_doughnuts Boolean. Plot the doughnuts that were created. Optional but highly
#'      highly recommended to verify the buffering does what you expect
#' @param object_n  The (row number of the) measurement point to plot. Only accepts single integer.
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


