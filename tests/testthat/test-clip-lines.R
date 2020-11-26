test_that("basic clipping", {
  x <- c(0, 0, 1, 1, 0)
  y <- c(0, 1, 1, 0, 0)
  id <- rep(1L, 5)

  # clip box entirely inside line
  out <- clip_lines_impl(x, y, id, .5, .5, .8, .8, 0)
  expect_identical(out$x, x)
  expect_identical(out$y, y)
  expect_identical(out$id, id)

  # clip box entirely outside line
  out <- clip_lines_impl(x, y, id, .5, .5, 1.2, 1.2, 0)
  expect_identical(out$x, numeric(0))
  expect_identical(out$y, numeric(0))
  expect_identical(out$id, integer(0))

  # clip top right corner
  out <- clip_lines_impl(x, y, id, 1, 1, 1, 1, 0)
  expect_equal(out$x, c(0.0, 0.0, 0.5, 1.0, 1.0, 0.0))
  expect_equal(out$y, c(0.0, 1.0, 1.0, 0.5, 0.0, 0.0))
  expect_identical(out$id, c(rep(1L, 3), rep(2L, 3)))

  # clip bottom left corner
  out <- clip_lines_impl(x, y, id, 0, 0, 1, 1, 0)
  expect_equal(out$x, c(0.0, 0.0, 1.0, 1.0, 0.5))
  expect_equal(out$y, c(0.5, 1.0, 1.0, 0.0, 0.0))
  expect_identical(out$id, rep(1L, 5))

  # clip right half
  out <- clip_lines_impl(x, y, id, 1, .5, 1, 2, 0)
  expect_equal(out$x, c(0.0, 0.0, 0.5, 0.5, 0.0))
  expect_equal(out$y, c(0, 1, 1, 0, 0))
  expect_identical(out$id, c(rep(1L, 3), 2L, 2L))

  # clip left half
  out <- clip_lines_impl(x, y, id, 0, .5, 1, 2, 0)
  expect_equal(out$x, c(0.5, 1.0, 1.0, 0.5))
  expect_equal(out$y, c(1, 1, 0, 0))
  expect_identical(out$id, rep(1L, 4))

  # clip in middle
  out <- clip_lines_impl(x, y, id, .5, .5, 2, .5, 0)
  expect_equal(out$x, c(0, 0, 0, 0, 1, 1, 1, 1, 0))
  expect_equal(out$y, c(0.00, 0.25, 0.75, 1.00, 1.00, 0.75, 0.25, 0.00, 0.00))
  expect_identical(out$id, c(rep(1L, 2), rep(2L, 4), rep(3L, 3)))

  out <- clip_lines_impl(x, y, id, .5, .5, .5, 2, 0)
  expect_equal(out$x, c(0.00, 0.00, 0.25, 0.75, 1.00, 1.00, 0.75, 0.25, 0.00))
  expect_equal(out$y, c(0, 1, 1, 1, 1, 0, 0, 0, 0))
  expect_identical(out$id, c(rep(1L, 3), rep(2L, 4), rep(3L, 2)))
})

test_that("clip multiple line segments", {
  x <- c(0, 0, 1, 1, 2, 2, 3, 3)
  y <- c(0, 1, 1, 0, 2, 3, 3, 2)
  id <- c(rep(1L, 4), rep(2L, 4))

  out <- clip_lines_impl(x, y, id, 0, 0, 1, 3, 0)
  expect_equal(out$x, c(0.5, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0))
  expect_equal(out$y, c(1, 1, 0, 2, 3, 3, 2))
  expect_identical(out$id, c(rep(1L, 3), rep(2L, 4)))

  out <- clip_lines_impl(x, y, id, 1.5, 1.5, 2, 4, 0)
  expect_equal(out$x, c(0.0, 0.0, 0.5, 2.5, 3.0, 3.0))
  expect_equal(out$y, c(0, 1, 1, 3, 3, 2))
  expect_identical(out$id, c(rep(1L, 3), rep(2L, 3)))
})

test_that("rotated clip box", {
  x <- c(0, 0, 1, 1, 0)
  y <- c(0, 1, 1, 0, 0)
  id <- rep(1L, 5)

  out <- clip_lines_impl(x, y, id, .5, .5, 5, sin(2*pi*45/360), 2*pi*45/360)
  expect_equal(out$x, c(0.0, 0.0, 0.5, 1.0, 1.0, 0.5))
  expect_equal(out$y, c(0.5, 1, 1, 0.5, 0, 0))
  expect_identical(out$id, c(rep(1L, 3), rep(2L, 3)))

  out <- clip_lines_impl(x, y, id, .5, .5, 5, sin(2*pi*45/360), -2*pi*45/360)
  expect_equal(out$x, c(0.0, 0.0, 0.5, 1.0, 1.0, 0.5, 0.0))
  expect_equal(out$y, c(0.0, 0.5, 1, 1, 0.5, 0, 0))
  expect_identical(out$id, c(rep(1L, 2), rep(2L, 3), rep(3L, 2)))

  out <- clip_lines_impl(x, y, id, .5, .5, sin(2*pi*45/360), 5, 2*pi*45/360)
  expect_equal(out$x, c(0.0, 0.0, 0.5, 1.0, 1.0, 0.5, 0.0))
  expect_equal(out$y, c(0.0, 0.5, 1, 1, 0.5, 0, 0))
  expect_identical(out$id, c(rep(1L, 2), rep(2L, 3), rep(3L, 2)))
})

test_that("singletons are carried over or clipped", {
  x <- c(0, 0, 1, 1, 2, 2, 3, 3)
  y <- c(0, 1, 1, 0, 2, 3, 3, 2)
  id <- c(1L, rep(2L, 3), rep(3L, 4))

  out <- clip_lines_impl(x, y, id, 10, 10, 1, 1, 0)
  expect_identical(out$x, x)
  expect_identical(out$y, y)
  expect_identical(out$id, id)

  out <- clip_lines_impl(x, y, id, 0, 0, .1, .1, 0)
  expect_identical(out$x, x[2:8])
  expect_identical(out$y, y[2:8])
  expect_identical(out$id, id[2:8] - 1L)

  x <- c(0, 0, 1, 1, 2, 2, 3, 3)
  y <- c(0, 1, 1, 0, 2, 3, 3, 2)
  id <- c(rep(1L, 3), 2L, rep(3L, 4))

  out <- clip_lines_impl(x, y, id, 10, 10, 1, 1, 0)
  expect_equal(out$x, x)
  expect_equal(out$y, y)
  expect_identical(out$id, id)

  out <- clip_lines_impl(x, y, id, 1, 0, .1, .1, 0)
  expect_identical(out$x, x[c(1:3, 5:8)])
  expect_identical(out$y, y[c(1:3, 5:8)])
  expect_identical(out$id, c(rep(1L, 3), rep(2L, 4)))

  x <- c(0, 0, 1, 1, 2, 2, 3, 3)
  y <- c(0, 1, 1, 0, 2, 3, 3, 2)
  id <- c(rep(1L, 3), rep(2L, 4), 3L)

  out <- clip_lines_impl(x, y, id, 10, 10, 1, 1, 0)
  expect_equal(out$x, x)
  expect_equal(out$y, y)
  expect_identical(out$id, id)

  out <- clip_lines_impl(x, y, id, 3, 2, .1, .1, 0)
  expect_identical(out$x, x[1:7])
  expect_identical(out$y, y[1:7])
  expect_identical(out$id, id[1:7])
})

test_that("empty or incorrect input", {
  out <- clip_lines_impl(numeric(0), numeric(0), integer(0), 3, 2, .1, .1, 0)
  expect_identical(out$x, numeric(0))
  expect_identical(out$y, numeric(0))
  expect_identical(out$id, integer(0))

  expect_error(
    clip_lines_impl(numeric(0), numeric(1), integer(0), 3, 2, .1, .1, 0),
    "must match"
  )
  expect_error(
    clip_lines_impl(numeric(0), numeric(0), integer(1), 3, 2, .1, .1, 0),
    "must match"
  )
  expect_error(
    clip_lines_impl(numeric(1), numeric(0), integer(0), 3, 2, .1, .1, 0),
    "must match"
  )
})
