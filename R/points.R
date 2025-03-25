#' Sum number of points within distance and polygon
#'
#' @inheritParams get_pts_buffer
#' @param dist The distance
#'
#' @returns list containing sum of points and area in buffer zone
#' @export
#'
sum_within_dist <- function(x, y, poly, dist){
  buff <- st_buffer(x, dist)
  y_in_buff <- st_contains(buff, y)
  sum_y_in_buff <- sapply(y_in_buff, length)
  suppressWarnings(area_in_buff <- st_area(st_intersection(buff, poly)))
  list(sum_y_in_buff,area_in_buff)
}


#' Get number of points within varying distances
#'
#' @param x The measurement points
#' @param y The points to be summed
#' @param poly The polygon to be kept within
#' @param mindist The minimum distance
#' @param maxdist The maximum distance
#' @param incdist The distance increment
#'
#' @returns list containing number of points and buffer area for every distance
#'   increment
#' @export
#'
get_pts_buffer <- function(x, y, poly, mindist, maxdist, incdist){
  dist <- seq(mindist,maxdist,incdist)
  sumdist <- lapply(dist, sum_within_dist, x = x, y = y, poly = poly)
  num_points <- sapply(sumdist, "[[", 1)
  buffer_area <- sapply(sumdist, "[[", 2)

  # convert from cumulative
  num_points <- t(apply(num_points, 1, function(x) diff(c(0,x))))
  buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,x))))
  buffer_area <- 1e-4*buffer_area #convert to hectares

  return(list(num_points = num_points,
              buffer_area = buffer_area))
}
