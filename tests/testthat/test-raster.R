# simulate data

set.seed(100230)
boundary <- sf::st_polygon(list(matrix(c(0,1000,1000,0,0,
                                         0,0,1000,1000,0),ncol=2)))
seedtraps <- data.frame(
  x = rep(seq(50,950,100), 10),
  y = rep(seq(50,950,100), each = 10)
) |>
  sf::st_as_sf(coords = c("x","y"), remove = FALSE,
               crs = "EPSG:32630")

dist_seq <- seq(10,200,length.out=10)

r <- terra::rast(terra::vect(boundary), nrow = 100, ncol = 100, vals = 1, crs = "EPSG:32630")
renv <- prioritizr::simulate_data(x = r, n = 1, scale = 0.8, intensity = 22, sd = 1)

growth_rate <- stats::rnorm(100)


# test run
test_that("raster data extraction works", {
  expect_error(lg_rast(x = seedtraps, y = renv, func = 2,
                        dist_seq = dist_seq, mindist = NULL, maxdist = NULL, incdist = NULL,
                        plot_doughnuts = FALSE),
               "'func' must be of class 'character' or 'function'")

  expect_message(rast_sum <- lg_rast(x = seedtraps, y = renv,
                                     dist_seq = c(10,20), mindist = NULL, maxdist = NULL, incdist = NULL,
                                     plot_doughnuts = TRUE, bg_layer = renv),
                 "plotting object indices")
  expect_named(rast_sum, c("buffer_values", "buffer_area", "value_parea"))

  rast_sum <- lg_rast(x = seedtraps,
                      y = renv,
                      func = "mean",
                      dist_seq = dist_seq, mindist = NULL, maxdist = NULL, incdist = NULL,
                      plot_doughnuts = FALSE)
  expect_named(rast_sum, c("buffer_values", "buffer_area", "value_parea"))
  expect_equal(dim(rast_sum$buffer_values), c(100,10))
  expect_equal(dim(rast_sum$buffer_area), c(100,10))
  expect_equal(dim(rast_sum$value_parea), c(100,10))

  lgr <- lg_assembly(response = growth_rate,
                     covariate = rast_sum$value_parea,
                     index_seq = dist_seq)
  expect_named(lgr, c("response", "covariate", "index"))
  expect_length(lgr$response, 100)
  expect_equal(dim(lgr$covariate), c(100,10))
  expect_equal(dim(lgr$index), c(100,10))

  pl <- lg_plot(lgr)
  expect_s3_class(pl, "patchwork")
  expect_length(pl, 3)
})

