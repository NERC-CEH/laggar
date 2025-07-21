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
  expect_named(res1, c("num_points", "buffer_area", "points_parea"))
  expect_true(all(is.na(res1[[2]])))
  expect_true(all(is.na(res1[[3]])))


  expect_no_error(res2 <- lg_points(test_x, test_y, seq(1,10,2),
                                    bound_poly = boundary))
  expect_named(res2, c("num_points", "buffer_area", "points_parea"))
  expect_false(any(is.na(res2[[2]])))
  expect_false(any(is.na(res2[[3]])))
})
