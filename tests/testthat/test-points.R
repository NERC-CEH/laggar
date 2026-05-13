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
                      dist_seq = dist_seq, bound_poly = NULL)
  expect_named(points, c("buffer_values", "buffer_area", "value_parea"))
  expect_equal(dim(points$buffer_values), c(100,20))
  expect_equal(points$buffer_area, NA)
  expect_equal(points$value_parea, NA)

  points_parea <- lg_points(seedtraps, adults,
                            dist_seq = dist_seq, bound_poly = boundary)
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


test_x <- sf::st_as_sf(data.frame(X = 1:40, Y = 40:1),
                        coords = c("X","Y"))
test_y <- sf::st_as_sf(data.frame(X = runif(100,1,40),
                                  Y = runif(100,1,40)),
                       coords = c("X","Y"))

boundary <- sf::st_polygon(list(matrix(c(0,45,45,0,0,
                                         0,0,45,45,0),ncol=2)))
boundary2 <- boundary + 5
comb <- sf::st_sfc(boundary, boundary2)

test_that("lg_points errors when expected", {
  expect_error(lg_points(test_x, boundary, seq(1,10,1)),
               "including only POINT geometries")
  expect_error(lg_points(boundary, test_y, seq(1,10,1)),
               "including only POINT geometries")
  expect_error(lg_points(test_x, test_y, seq(1,10,1),
                         bound_poly = comb),
               "sf object comprising a single polygon")
  expect_error(lg_points(test_x, test_y, seq(1,10,1),
                         bound_poly = test_x),
               "sf object comprising a single polygon")
})

test_that("lg_points runs", {
  expect_no_error(res1 <- lg_points(test_x, test_y, seq(1,10,1)))
  expect_named(res1, c("buffer_values", "buffer_area", "value_parea"))
  expect_true(all(is.na(res1[[2]])))
  expect_true(all(is.na(res1[[3]])))


  expect_no_error(res2 <- lg_points(test_x, test_y, seq(1,10,2),
                                    bound_poly = boundary))
  expect_named(res2, c("buffer_values", "buffer_area", "value_parea"))
  expect_false(any(is.na(res2[[2]])))
  expect_false(any(is.na(res2[[3]])))
})

