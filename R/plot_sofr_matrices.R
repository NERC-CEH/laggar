#' Visualize scalar-on-function regression data
#'
#' Plots the response, covariate matrix and index matrix for a given
#' SoFR term in a model.
#'
#' @param data a `data.frame`
#' @param response string giving the name of the column with the response in it
#' @param covariate string giving the name of the matrix-column with covariate values in it
#' @param index string giving the name of the matrix-column with index values in it
#' @return a `ggplot2` plot
#'
#' @author David L Miller
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate "%>%"
#' @importFrom rlang "!!"
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_tile aes_string scale_fill_viridis_c scale_x_continuous expansion scale_y_continuous labs theme_minimal theme element_blank
#' @importFrom patchwork wrap_plots
plot_sofr_matrices <- function(data, response, covariate, index){

  # drop geometries
  data <- st_drop_geometry(data)

  # re-arrange the data to be ggplot2-compatible
  covm <- as.data.frame(data[[covariate]])
  indm <- as.data.frame(data[[index]])

  data$.sample_id <- 1:nrow(data)

  colnames(covm) <- paste0(covariate,1:ncol(covm))
  colnames(indm) <- paste0(index,1:ncol(indm))

  covdata <- cbind(data, covm) %>%
    select(!!colnames(covm),  !!response, ".sample_id") %>%
    pivot_longer(!c(!!response, ".sample_id"),
                 values_to=covariate, names_to="mo") %>%
    mutate(mo = as.numeric(sub(covariate, "", mo)))
  inddata <- cbind(data, indm) %>%
    select(!!colnames(indm),  !!response, ".sample_id") %>%
    pivot_longer(!c(!!response, ".sample_id"),
                 values_to=index, names_to="mo") %>%
    mutate(mo = as.numeric(sub(index, "", mo)))

  data$.x <- 1

  pl <- list(ggplot(data) +
              geom_tile(aes_string(x=".x", y=".sample_id",
                                   fill=response)) +
              scale_fill_viridis_c() +
              scale_x_continuous(expand=expansion(0)) +
              scale_y_continuous(expand=expansion(0)) +
              labs(x=response, y="Sample") +
              theme_minimal() +
              theme(legend.position="bottom",
                    axis.text.x=element_blank()),
             ggplot(covdata) +
              geom_tile(aes_string(x="mo", y=".sample_id",
                                   fill=covariate)) +
              scale_fill_viridis_c() +
              scale_x_continuous(expand=expansion(0)) +
              scale_y_continuous(expand=expansion(0)) +
              labs(x=index, y="") +
              theme_minimal() +
              theme(legend.position="bottom",
                    axis.text.y=element_blank()),
             ggplot(inddata) +
              geom_tile(aes_string(x="mo", y=".sample_id",
                                   fill=index)) +
              scale_fill_viridis_c() +
              scale_x_continuous(expand=expansion(0)) +
              scale_y_continuous(expand=expansion(0)) +
              labs(x=index, y="") +
              theme_minimal() +
              theme(legend.position="bottom",
                    axis.text.y=element_blank()))

  wrap_plots(pl, ncol=3, widths=c(0.1, 0.5, 0.5))
}
