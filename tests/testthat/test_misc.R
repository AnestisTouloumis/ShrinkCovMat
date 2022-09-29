test_that("centered test", {
  centered <- as.logical("some character")
  expect_error(if (centered != TRUE && centered != FALSE) {
    stop("'centered' must be either 'TRUE' or 'FALSE'")
  })
  centered <- TRUE
  expect_false(centered != TRUE)
  expect_true(centered != FALSE)
  centered <- FALSE
  expect_true(centered != TRUE)
  expect_false(centered != FALSE)
  n <- 5
  expect_false(n < 2)
  expect_false(n < 4)
  n <- 1
  expect_error(
    if (n < 2) stop("The number of columns should be greater than 1")
  )
  n <- sample(1:3, 1)
  expect_error(
    if (n < 4) stop("The number of columns should be greater than 3")
  )
})
