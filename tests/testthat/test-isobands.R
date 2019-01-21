context("test-isobands")

test_that("elementary polygons get merged", {
  # two connected polygons get merged
  z <- matrix(c(0, 0, 1,
                1, 1, 1), ncol = 3, nrow = 2, byrow = TRUE)
  out <- isobands(x = 1:3, y = 2:1, z, levels_low = 0.5, levels_high = 1.5)
  expect_setequal(10*out[[1]]$x + out[[1]]$y,
                  10*c(3.0, 2.0, 1.0, 1.0, 2.0, 2.5, 3.0) +
                    c(1.0, 1.0, 1.0, 1.5, 1.5, 2.0, 2.0))
  expect_equal(out[[1]]$id, rep(1, 7))
  #
  # two unconnected polygons don't get merged
  z <- matrix(c(1, 2, 1,
                1, 2, 2), ncol = 3, nrow = 2, byrow = TRUE)
  out <- isobands(x = 1:3, y = 2:1, z, levels_low = 0.5, levels_high = 1.5)
  expect_setequal(10*out[[1]]$x + out[[1]]$y,
                  10*c(3.0, 2.5, 3.0, 1.0, 1.5, 1.5, 1.0) +
                    c(1.5, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0))
  expect_setequal(out[[1]]$id, c(1:2))
  expect_equal(length(out[[1]]$id), 7)
  #
  # # two separate lines get merged in second row
  # z <- matrix(c(0, 1, 0,
  #               0, 1, 0,
  #               0, 0, 0), ncol = 3, nrow = 3, byrow = TRUE)
  # out <- isolines(x = 1:3, y = 3:1, z, levels = 0.5)
  # expect_setequal(out[[1]]$x, c(2.5, 2.5, 2.0, 1.5, 1.5))
  # expect_setequal(out[[1]]$y, c(3.0, 2.0, 1.5, 2.0, 3.0))
  # expect_equal(out[[1]]$id, rep(1, 5))
  #
  # # circle gets closed
  # z <- matrix(c(0, 0, 0,
  #               0, 1, 0,
  #               0, 0, 0), ncol = 3, nrow = 3, byrow = TRUE)
  # out <- isolines(x = 1:3, y = 3:1, z, levels = 0.5)
  # # circle is closed
  # expect_equal(out[[1]]$x[1], out[[1]]$x[5])
  # expect_equal(out[[1]]$y[1], out[[1]]$y[5])
  # # coords are correct
  # expect_setequal(out[[1]]$x, c(2.5, 2.0, 1.5, 2.0, 2.5))
  # expect_setequal(out[[1]]$y, c(2.0, 2.5, 2.0, 1.5, 2.0))
  # expect_equal(out[[1]]$id, rep(1, 5))
})
