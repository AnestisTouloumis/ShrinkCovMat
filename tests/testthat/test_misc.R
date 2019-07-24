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
})
