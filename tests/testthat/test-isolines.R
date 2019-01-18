context("test-isolines")

test_that("line segments get merged", {
  # two connected line segments get merged
  z <- matrix(c(0, 0, 1,
                1, 1, 1), ncol = 3, nrow = 2, byrow = TRUE)
  out <- isolines(x = 1:3, y = 2:1, z, levels = 0.5)
  expect_equal(out[[1]]$x, c(1, 2, 2.5))
  expect_equal(out[[1]]$y, c(1.5, 1.5, 2.0))
  expect_equal(out[[1]]$id, rep(1, 3))

  # two unconnected line segments don't get merged
  z <- matrix(c(0, 1, 0,
                0, 1, 1), ncol = 3, nrow = 2, byrow = TRUE)
  out <- isolines(x = 1:3, y = 2:1, z, levels = 0.5)
  expect_equal(out[[1]]$x, c(2.5, 3.0, 1.5, 1.5))
  expect_equal(out[[1]]$y, c(2.0, 1.5, 2.0, 1.0))
  expect_equal(out[[1]]$id, rep(1:2, each = 2))

  # two separate lines get merged in second row
  z <- matrix(c(0, 1, 0,
                0, 1, 0,
                0, 0, 0), ncol = 3, nrow = 3, byrow = TRUE)
  out <- isolines(x = 1:3, y = 3:1, z, levels = 0.5)
  expect_equal(out[[1]]$x, c(2.5, 2.5, 2.0, 1.5, 1.5))
  expect_equal(out[[1]]$y, c(3.0, 2.0, 1.5, 2.0, 3.0))
  expect_equal(out[[1]]$id, rep(1, 5))

  # circle gets closed
  z <- matrix(c(0, 0, 0,
                0, 1, 0,
                0, 0, 0), ncol = 3, nrow = 3, byrow = TRUE)
  out <- isolines(x = 1:3, y = 3:1, z, levels = 0.5)
  expect_equal(out[[1]]$x, c(2.5, 2.0, 1.5, 2.0, 2.5))
  expect_equal(out[[1]]$y, c(2.0, 2.5, 2.0, 1.5, 2.0))
  expect_equal(out[[1]]$id, rep(1, 5))

})
