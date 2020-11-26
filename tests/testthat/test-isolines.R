test_that("line segments get merged", {
  # two connected line segments get merged
  z <- matrix(c(0, 0, 1,
                1, 1, 1), ncol = 3, nrow = 2, byrow = TRUE)
  out <- isolines(x = 1:3, y = 2:1, z, levels = 0.5)
  expect_setequal(10000*out[[1]]$x + out[[1]]$y,
                  10000*c(1, 2, 2.5) + c(1.5, 1.5, 2.0))
  expect_equal(out[[1]]$id, rep(1, 3))

  # two unconnected line segments don't get merged
  z <- matrix(c(0, 1, 0,
                0, 1, 1), ncol = 3, nrow = 2, byrow = TRUE)
  out <- isolines(x = 1:3, y = 2:1, z, levels = 0.5)
  expect_setequal(10000*out[[1]]$x + out[[1]]$y,
                  10000*c(2.5, 3.0, 1.5, 1.5) + c(2.0, 1.5, 2.0, 1.0))
  expect_setequal(out[[1]]$id, c(1:2))
  expect_equal(length(out[[1]]$id), 4)

  # two separate lines get merged in second row
  z <- matrix(c(0, 1, 0,
                0, 1, 0,
                0, 0, 0), ncol = 3, nrow = 3, byrow = TRUE)
  out <- isolines(x = 1:3, y = 3:1, z, levels = 0.5)
  expect_setequal(10000*out[[1]]$x + out[[1]]$y,
                  10000*c(2.5, 2.5, 2.0, 1.5, 1.5) + c(3.0, 2.0, 1.5, 2.0, 3.0))
  expect_equal(out[[1]]$id, rep(1, 5))

  # circle gets closed
  z <- matrix(c(0, 0, 0,
                0, 1, 0,
                0, 0, 0), ncol = 3, nrow = 3, byrow = TRUE)
  out <- isolines(x = 1:3, y = 3:1, z, levels = 0.5)
  # circle is closed
  expect_equal(out[[1]]$x[1], out[[1]]$x[5])
  expect_equal(out[[1]]$y[1], out[[1]]$y[5])
  # coords are correct
  expect_setequal(10000*out[[1]]$x + out[[1]]$y,
                  10000*c(2.5, 2.0, 1.5, 2.0, 2.5) + c(2.0, 2.5, 2.0, 1.5, 2.0))
  expect_equal(out[[1]]$id, rep(1, 5))
})

test_that("NAs are handled correctly", {
  z <- matrix(c(NA, 0, 0,
                 0, 1, 1,
                 0, 1, 1), ncol = 3, nrow = 3, byrow = TRUE)
  out <- isolines(x = 1:3, y = 3:1, z, levels = 0.5)
  expect_setequal(10000*out[[1]]$x + out[[1]]$y,
                  10000*c(1.5, 1.5, 2.0, 3.0) +
                    c(2.0, 1.0, 2.5, 2.5))
  expect_setequal(out[[1]]$id, c(1:2))
  expect_equal(length(out[[1]]$id), 4)
})


test_that("All elementary segments are calculated correctly", {
  # a matrix that requires all elementary segments for isolines
  z <- matrix(c(0, 0, 0, 1, 1, 0, 1, 1,
                0, 0, 0, 1, 1, 0, 1, 1,
                0, 1, 1, 0, 1, 1, 0, 0,
                1, 1, 0, 0, 0, 1, 1, 0,
                1, 0, 1, 1, 0, 0, 0, 1
  ), ncol = 8, nrow = 5, byrow = TRUE)
  out <- isolines(x = 1:8, y = 5:1, z, levels = 0.5)

  expect_setequal(
    10000*out[[1]]$x + out[[1]]$y,
    10000*c(7.5, 7.0, 6.0, 5.5, 5.0, 4.5, 4.0, 3.5,
            3.0, 2.5, 3.0, 4.0, 4.5, 2.5, 2.0, 1.5,
            8.0, 7.0, 6.5, 7.0, 7.5, 8.0, 6.5, 6.5,
            6.0, 5.5, 5.5, 3.5, 3.5, 3.0, 2.0, 1.5,
            1.0) +
          c(1.0, 1.5, 1.5, 2.0, 2.5, 3.0, 3.5, 3.0,
            2.5, 2.0, 1.5, 1.5, 1.0, 1.0, 1.5, 1.0,
            3.5, 3.5, 3.0, 2.5, 2.0, 1.5, 5.0, 4.0,
            3.5, 4.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,
            2.5)
  )
  expect_setequal(out[[1]]$id, c(1:5))
  expect_equal(length(out[[1]]$id), 33)
})


test_that("Saddles", {
  # a matrix that contains all saddles (there are only two)
  z <- matrix(c(0, 1, 0,
                1, 0, 1),
              ncol = 3, nrow = 2, byrow = TRUE)

  out <- isolines(x = 1:3, y = 2:1, z, levels = 0.5)

  expect_setequal(10000*out[[1]]$x + out[[1]]$y,
                  10000*c(2.5, 3.0, 2.5, 2.0, 1.5, 1.5, 1.0) +
                        c(2.0, 1.5, 1.0, 1.5, 1.0, 2.0, 1.5))
  expect_setequal(out[[1]]$id, c(1:3))
  expect_equal(length(out[[1]]$id), 7)

  out <- isolines(x = 1:3, y = 2:1, z, levels = 0.6)

  expect_setequal(10000*out[[1]]$x + out[[1]]$y,
                  10000*c(1.6, 2.0, 2.4, 3.0, 2.6, 1.0, 1.4) +
                        c(2.0, 1.6, 2.0, 1.4, 1.0, 1.4, 1.0))
  expect_setequal(out[[1]]$id, c(1:3))
  expect_equal(length(out[[1]]$id), 7)
})
