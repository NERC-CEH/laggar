test_that("lg_assembly errors correctly", {
  expect_error(lg_assembly(c("abc","def"), matrix(1:20,nrow=4),
                           index_seq = 1:5),
               "response must be a numeric")

  expect_error(lg_assembly(1:5, matrix(1:20,nrow=4), index_seq=1:5),
               "number of rows")

  expect_error(lg_assembly(1:4, as.data.frame(matrix(1:20,nrow=4)),
                           index_seq=1:5),
               "covariate must be a matrix")
})

test_that("lg_assembly handles NAs correctly", {
  expect_error(lg_assembly(1:4, matrix(c(1:16,rep(NA,4)), nrow=4),
                           index_seq = 1:6),
               "has NAs present in every row")

  expect_error(lg_assembly(as.numeric(rep(NA, 4)), matrix(1:20, nrow=4),
                           index_seq = 1:6),
               "Response is entirely NAs")

  expect_message(a1 <- lg_assembly(1:4, matrix(c(1:19,NA), nrow=4),
                                   index_seq = seq(1,5,1)),
                 "Removing rows where covariate matrix has NAs")
  expect_named(a1, c("response", "covariate", "index"))
  expect_length(a1$response, 3)
  expect_equal(dim(a1$covariate),c(3,5))
  expect_equal(dim(a1$index),c(3,5))


  expect_message(a2 <- lg_assembly(c(1,NA,NA,4), matrix(1:20, nrow=4),
                                   index_seq = seq(1,5,1)),
                 "Removing rows where response is NA")
  expect_named(a2, c("response", "covariate", "index"))
  expect_length(a2$response, 2)
  expect_equal(dim(a2$covariate),c(2,5))
  expect_equal(dim(a2$index),c(2,5))
})

# Helper function tests ####
test_that("index checking errors as expected", {
  expect_error(.index_check(0,200,NULL,NULL),
               "must be specified")
  expect_error(.index_check(-20,5,6,NULL),
               "must be greater than 0")
  expect_error(.index_check(50,10,5,NULL),
               "must be greater than")
  expect_error(.index_check(NULL,NULL,NULL,
                            seq(-5,10,1)),
               "must be greater than 0")

  expect_error(.index_check("a","b","c",NULL),
               "must be numeric")
  expect_error(.index_check("a","b","c","d"),
               "must be numeric")
  expect_error(.index_check(NULL, NULL, NULL,1),
               "must be of length greater than 1")

})

test_that("index checking works as expected", {
  expect_no_error(.index_check(NULL,NULL,NULL,
                               seq(10,200,10)))
  expect_length(.index_check(1,5,1,NULL), 5)
  expect_length(.index_check(11,213,9,NULL),23)

})


test_sf <- sf::st_as_sf(data.frame(X = 1:40, Y = 40:1),
                        coords = c("X","Y"))
# test_sp <- sf::as_Spatial(test_sf) # removing as don't want dependency on sp

test_that("sf conversion errors as expected", {
  expect_error(suppressMessages(.sf_conversion("x_sp")),
               "no applicable method")
})
test_that("sf conversion works as expected", {
  # expect_message(test2 <- .sf_conversion(test_sp),
  #                "Converting object")
  # expect_s3_class(test2, "sf")

  expect_no_message(test3 <- .sf_conversion(test_sf))
  expect_s3_class(test3, "sf")
})
