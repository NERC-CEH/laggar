
test_that("lg_assembly errors when expected", {
  expect_error(lg_assembly(response = LETTERS[1:10],
                           index_seq = seq(1,2,1), covariate = matrix(ncol=2)),
               "response must be a numeric or integer vector")
  expect_error(lg_assembly(response = rnorm(11),
                           index_seq = seq(1,3,1),
                           covariate = as.data.frame(matrix(rnorm(11*3),nrow=11))),
               "covariate must be a matrix")
  expect_error(lg_assembly(response = rnorm(11),
                           index_seq = seq(1,11,1),
                           covariate = matrix(rnorm(11*3),nrow=3)),
               "number of rows of covariate must be equal to length of response")
  expect_error(lg_assembly(response = rnorm(11),
                           index_seq = seq(1,4,1),
                           covariate = cbind(matrix(rnorm(11*3),nrow=11),NA)),
               "covariate matrix has NAs present in every row")
  expect_error(lg_assembly(response = as.numeric(rep(NA, 11)),
                           index_seq = seq(1,3,1),
                           covariate = matrix(rnorm(11*3),nrow=11)),
               "Response is entirely NAs")
})

test_that("lg_assembly handles NAs okay", {
  expect_warning(
    suppressMessages(
      lg1 <- lg_assembly(response = seq(1,12),
                         index_seq = seq(1,10,1),
                         covariate = rbind(matrix(rnorm(11*9),nrow=11),NA))),
    "Index sequence must be of the same length as the number of columns"
  )
  expect_equal(lg1$response, seq(1,11))
  expect_equal(dim(lg1$index), c(11,9))
  expect_equal(dim(lg1$covariate), c(11,9))
  expect_all_true(apply(lg1$index, 2, function(x) length(unique(x))==1))
  expect_all_true(apply(lg1$index, 1, function(x) all.equal(x, 1:9)))


  expect_message(
    lg1 <- lg_assembly(response = seq(1,12),
                       index_seq = seq(1,9,1),
                       covariate = rbind(matrix(rnorm(11*9),nrow=11),NA)),
    "Removing rows where covariate matrix has NAs")
  expect_named(lg1, c("response", "covariate", "index"))
  expect_equal(lg1$response, seq(1,11))
  expect_equal(dim(lg1$index), c(11,9))
  expect_equal(dim(lg1$covariate), c(11,9))


  expect_message(
    lg1 <- lg_assembly(response = c(seq(1,10),NA),
                       index_seq = seq(1,9,1),
                       covariate = matrix(rnorm(11*9),nrow=11)),
    "Removing rows where response is NA")
  expect_named(lg1, c("response", "covariate", "index"))
  expect_equal(lg1$response, seq(1,10))
  expect_equal(dim(lg1$index), c(10,9))
  expect_equal(dim(lg1$covariate), c(10,9))

})

# helper function checks
test_that("index_check errors when expected", {
  expect_error(.index_check(3,4,NULL,index_seq = NULL),
               "Either index_seq or minindex, maxindex and incindex must be specified")
  expect_error(.index_check("a","b","c",NULL),
               "numeric")
  expect_error(.index_check(c(1,2),c(1,3),c(1,2),NULL),
               "must be of length")
  expect_error(.index_check(-5,3,5,NULL),
               "must be greater than 0")
  expect_error(.index_check(5,4,1,NULL),
               "must be greater than")
  expect_error(.index_check(NULL,NULL,NULL,index_seq = "aw"),
               "numeric")
  expect_error(.index_check(NULL,NULL,NULL,index_seq = 1),
               "must be of length")
  expect_error(.index_check(NULL,NULL,NULL,index_seq = c(-5,2,3)),
               "must be greater than")
})
