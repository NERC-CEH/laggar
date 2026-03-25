# simulate data
set.seed(100230)
boundary <- sf::st_polygon(list(matrix(c(0,1000,1000,0,0,
                                         0,0,1000,1000,0),ncol=2)))
adults <- sf::st_sample(boundary, size = 300)
seedtraps <- data.frame(
  x = rep(seq(50,950,100), 10),
  y = rep(seq(50,950,100), each = 10)
) |>
  sf::st_as_sf(coords = c("x","y"), remove = FALSE)

nseeds <- stats::rpois(100, 1)
dist_seq <- seq(10,200,10)

# test run
test_that("points data extraction works", {
  points <- lg_points(seedtraps, adults,
                      dist_seq = dist_seq, poly = NULL)
  expect_named(points, c("buffer_values", "buffer_area", "value_parea"))
  expect_equal(dim(points$buffer_values), c(100,20))
  expect_equal(points$buffer_area, NA)
  expect_equal(points$value_parea, NA)

  points_parea <- lg_points(seedtraps, adults,
                            dist_seq = dist_seq, poly = boundary)
  expect_named(points_parea, c("buffer_values", "buffer_area", "value_parea"))
  expect_equal(dim(points_parea$buffer_values), c(100,20))
  expect_equal(dim(points_parea$buffer_area), c(100,20))
  expect_equal(dim(points_parea$value_parea), c(100,20))

  lgr <- lg_assembly(response = nseeds,
                     covariate = points_parea$value_parea,
                     index_seq = dist_seq)
  expect_named(lgr, c("response", "covariate", "index"))
  expect_length(lgr$response, 100)
  expect_equal(dim(lgr$covariate), c(100,20))
  expect_equal(dim(lgr$index), c(100,20))

  pl <- lg_plot(lgr)
  expect_s3_class(pl, "patchwork")
  expect_length(pl, 3)
})

