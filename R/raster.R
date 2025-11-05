

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
                             nobjects, # number of objects that were originally buffered
                             object_n = 1, # the number of object that you want to check (i.e. which buffered point do you want plotted). Can be NULL if want to skip single object check, not recommended.
                             bg_layer = NULL) # background raster layer for more plotting context
{

  if(nbuffers>9)
    stop("Maximum number of buffers for plotting is 9.")
  if(object_n > nobjects)
    stop("The object chosen for plotting must be in range 0:nobjects")


  # create an index for plotting
  index_for_plotting <-
    (0:(nbuffers-1) * nobjects) + object_n
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
    ggplot() +
      {if(!is.null(bg_layer)) geom_tile(data = bg_layer, aes(x, y, fill = val))} +
      scale_fill_viridis_c(option = "D", name = NULL) +
      ggnewscale::new_scale_fill() +
      geom_sf(data = doughnut_out[index_for_plotting[x],],
              aes(fill =  buffer_name), alpha = 0.7) +

      scale_fill_manual(values = colpal, name = NULL, guide = "none") +
      ggtitle(paste("Buffer distance", unique(doughnut_out[index_for_plotting[x],]$buffer_name))) +
      coord_sf(xlim = st_bbox(doughnut_out[max(index_for_plotting),])[c(1, 3)],
               ylim = st_bbox(doughnut_out[max(index_for_plotting),])[c(2, 4)]) +
      theme_bw()
  )

  patchwork::wrap_plots(plts, guides = "collect")

}


