
test_resp <- rnorm(10)
test_cov <- matrix(rnorm(10*5),nrow=10)
test_ind <- matrix(rep(1:5,each=10),nrow=10)

test_that("lg_plot.list works", {
  pl <- lg_plot(list(response = test_resp,
                     covariate = test_cov,
                     index = test_ind),
                widths = c(1,1,1))
  expect_s3_class(pl, "patchwork")
  expect_length(pl, 2)

  pl2 <- lg_plot(list(resp = test_resp,
                      cov = test_cov,
                      ind = test_ind),
                 response = "resp",
                 covariate = "cov",
                 index = "ind",
                 plot_index = TRUE)
  expect_s3_class(pl2, "patchwork")
  expect_length(pl2, 3)

  expect_error(lg_plot(list(response = test_resp,
                            covariate = test_cov,
                            index = test_ind),
                       response = "resp"),
               "not a named part")
})

test_df <- data.frame(response = test_resp)
test_df$covariate <- test_cov
test_df$index <- test_ind
test_that("lg_plot.data.frame works", {
  pl <- lg_plot(test_df,
                widths = c(1,1,1),
                plot_index = TRUE)
  expect_s3_class(pl, "patchwork")
  expect_length(pl, 3)

  expect_error(lg_plot(test_df,
                       response = "resp",
                       covariate = "cov",
                       index = "ind"),
               "not a named part")
})

test_sf <- sf::st_as_sf(dplyr::mutate(test_df, X = runif(10), Y = runif(10)),
                        coords = c("X","Y"))
test_that("lg_plot.sf works", {
  pl <- lg_plot(test_sf,
                widths = c(1,1,1),
                index = "dontfail")
  expect_s3_class(pl, "patchwork")
  expect_length(pl, 2)


  expect_error(lg_plot(test_sf,
                       index = "fail", plot_index = TRUE),
               "not a named part")


})
