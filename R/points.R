#' Sum number of points within distance and polygon
#'
#' @inheritParams lg_points
#' @param dist The distance
#'
#' @returns list containing sum of points and area in buffer zone
#' @export
#'
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


#' Get number of points within varying distances
#'
#' @param x The measurement points
#' @param y The points to be summed
#' @param dist_seq The sequence of distance values
#' @param poly Optional polygon to be kept within
#' @param mindist The minimum distance, ignored if dist_seq is specified
#' @param maxdist The maximum distance, ignored if dist_seq is specified
#' @param incdist The distance increment, ignored if dist_seq is specified
#'
#' @returns list containing matrices of the number of points
#'   (\code{buffer_values}), buffer area (\code{buffer_area}) and average number of
#'   points per unit area  (\code{value_parea}) for every distance increment.
#'   All matrices will have number of rows equal to the number of measurement
#'   points unless \code{poly = NULL} in which case \code{buffer_area} and
#'   \code{value_parea} will be returned as \code{NA}.
#' @export
#'
lg_points <- function(x, y, dist_seq = NULL, poly = NULL,
                      mindist = NULL, maxdist = NULL, incdist = NULL){
  dist <- .index_check(minindex = mindist, maxindex = maxdist,
                       incindex = incdist, index_seq = dist_seq)
  sumdist <- lapply(dist, sum_within_dist, x = x, y = y, poly = poly)
  num_points <- sapply(sumdist, "[[", 1)

  # convert from cumulative
  num_points <- t(apply(num_points, 1, function(x) diff(c(0,x))))
  if(!is.null(poly)){
    buffer_area <- sapply(sumdist, "[[", 2)
    buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,x))))
    points_parea <- num_points/buffer_area
  } else{
    buffer_area <- NA
    points_parea <- NA
  }

  return(list(buffer_values = num_points,
              buffer_area = buffer_area,
              value_parea = points_parea))
}

