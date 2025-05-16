#' Visualise scalar-on-function regression data
#'
#' Plots the response, covariate matrix and index matrix for a given SoFR term
#'
#' @param data Either a data.frame, sf object, or a list
#' @param response string giving the name of the column with the response in it
#' @param covariate string giving the name of the matrix-column with covariate values in it
#' @param index string giving the name of the matrix-column with index values in it
#' @param widths the widths of the different plots
#' @param ... Other arguments passed to [patchwork::wrap_plots()]
#'
#'
#' @return a `patchwork` combination of `ggplot2` objects
#'
#' @export
lg_plot <- function(...) UseMethod("lg_plot")

#' @describeIn lg_plot sf interface
#' @export
lg_plot.sf <- function(data, ...){
  # drop geometries
  data <- sf::st_drop_geometry(data)
  lg_plot.data.frame(data, ...)
}

#' @describeIn lg_plot data.frame interface
#' @export
#' @importFrom rlang .data
lg_plot.data.frame <- function(data, response, covariate, index,
                               widths = c(0.1, 0.5, 0.5), ...){

  # re-arrange the data to be ggplot2-compatible
  covm <- as.data.frame(data[[covariate]])
  indm <- as.data.frame(data[[index]])

  data$.sample_id <- 1:nrow(data)

  colnames(covm) <- paste0(covariate,1:ncol(covm))
  colnames(indm) <- paste0(index,1:ncol(indm))

  covdata <- cbind(data, covm) |>
    dplyr::select(.data[[colnames(covm)]],  .data[[response]], ".sample_id") |>
    tidyr::pivot_longer(!c(.data[[response]], ".sample_id"),
                 values_to=covariate, names_to="mo") |>
    dplyr::mutate(mo = as.numeric(sub(covariate, "", .data$mo)))
  inddata <- cbind(data, indm) |>
    dplyr::select(.data[[colnames(indm)]],  .data[[response]], ".sample_id") |>
    tidyr::pivot_longer(!c(.data[[response]], ".sample_id"),
                 values_to=index, names_to="mo") |>
    dplyr::mutate(mo = as.numeric(sub(index, "", .data$mo)))

  data$.x <- 1


  .make_plot(data = data, response = response,
             covdata = covdata, covariate = covariate,
             inddata = inddata, index = index, widths = widths, ...)
}

#' @describeIn lg_plot list interface
#' @export
lg_plot.list <- function(data, response = "response",
                         covariate = "covariate",
                         index = "index",
                         widths = c(0.1, 0.5, 0.5), ...){

  # re-arrange the data to be ggplot2-compatible
  covdata <- as.data.frame(data[[covariate]])
  inddata <- as.data.frame(data[[index]])
  respdf <- data.frame(resp = data[[response]])

  colnames(covdata) <- paste0(covariate,1:ncol(covdata))
  colnames(inddata) <- paste0(index,1:ncol(inddata))

  covdata$.sample_id <- 1:nrow(respdf)
  inddata$.sample_id <- 1:nrow(respdf)
  respdf$.sample_id <- 1:nrow(respdf)

  covdata <- tidyr::pivot_longer(covdata, !.data$.sample_id,
                                 values_to = covariate, names_to = "mo")
  covdata$mo <- as.numeric(sub(covariate, "", covdata$mo))
  inddata <- tidyr::pivot_longer(inddata, !.data$.sample_id,
                                 values_to = index, names_to = "mo")
  inddata$mo <- as.numeric(sub(index, "", inddata$mo))

  respdf$.x <- 1
  .make_plot(respdata = respdf, response = "resp",
             covdata = covdata, covariate = covariate,
             inddata = inddata, index = index,
             widths = widths, ...)
}

.make_plot <- function(respdata, response,
                       covdata, covariate,
                       inddata, index, widths, ...){
  pl <- list(ggplot2::ggplot(respdata) +
               ggplot2::geom_tile(ggplot2::aes_string(x=".x", y=".sample_id",
                                                      fill=response)) +
               ggplot2::scale_fill_viridis_c() +
               ggplot2::scale_x_continuous(expand=ggplot2::expansion(0)) +
               ggplot2::scale_y_continuous(expand=ggplot2::expansion(0)) +
               ggplot2::labs(x=response, y="Sample") +
               ggplot2::theme_minimal() +
               ggplot2::theme(legend.position="bottom",
                              axis.text.x=ggplot2::element_blank()),
             ggplot2::ggplot(covdata) +
               ggplot2::geom_tile(ggplot2::aes_string(x="mo", y=".sample_id",
                                                      fill=covariate)) +
               ggplot2::scale_fill_viridis_c() +
               ggplot2::scale_x_continuous(expand=ggplot2::expansion(0)) +
               ggplot2::scale_y_continuous(expand=ggplot2::expansion(0)) +
               ggplot2::labs(x=index, y="") +
               ggplot2::theme_minimal() +
               ggplot2::theme(legend.position="bottom",
                              axis.text.y=ggplot2::element_blank()),
             ggplot2::ggplot(inddata) +
               ggplot2::geom_tile(ggplot2::aes_string(x="mo", y=".sample_id",
                                                      fill=index)) +
               ggplot2::scale_fill_viridis_c() +
               ggplot2::scale_x_continuous(expand=ggplot2::expansion(0)) +
               ggplot2::scale_y_continuous(expand=ggplot2::expansion(0)) +
               ggplot2::labs(x=index, y="") +
               ggplot2::theme_minimal() +
               ggplot2::theme(legend.position="bottom",
                              axis.text.y=ggplot2::element_blank()))

  patchwork::wrap_plots(pl, ncol=3, widths=widths, ...)
}
