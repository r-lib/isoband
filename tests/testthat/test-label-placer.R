test_that("minmax label placer", {
  lines <- list(
    "1" = list(x = c(.5, .75, .5, .25), y = c(.25, .5, .75, .5), id = rep(1, 4)),
    "2" = list(x = c(.5, 1, .5, 0), y = c(0, .5, 1, .5), id = rep(1, 4))
  )

  labels_data <- data.frame(
    index = 1:2,
    break_index = c(3, 7),
    break_id = c("1", "2"),
    label = c("a", "b"),
    stringsAsFactors = FALSE
  )

  lp <- label_placer_minmax(n = 0)
  out <- lp(lines, labels_data)
  expect_equal(out$index, rep(1:2, each = 2))
  expect_equal(out$break_index, rep(c(3, 7), each = 2))
  expect_equal(out$break_id, rep(c("1", "2"), each = 2))
  expect_equal(out$label, rep(c("a", "b"), each = 2))
  expect_equal(out$x, rep(0.5, 4))
  expect_equal(out$y, c(0.75, 0.25, 1, 0))
  expect_equal(out$theta, rep(0, 4))

  lp <- label_placer_minmax(placement = "rl", rot_adjuster = angle_fixed(1.5), n = 0)
  out <- lp(lines, labels_data)
  expect_equal(out$x, c(0.25, 0.75, 0, 1))
  expect_equal(out$y, rep(0.5, 4))
  expect_equal(out$theta, rep(1.5, 4))

  lp <- label_placer_minmax(placement = NULL)
  out <- lp(lines, labels_data)
  expect_equal(nrow(out), 0)
})



test_that("angle adjustments", {
  theta <- c(-3, -2, -1, 0, 1, 2, 3)

  expect_equal(
    angle_halfcircle_bottom()(theta),
    c(-3 + pi, -2 + pi, -1, 0, 1, 2 - pi, 3 - pi)
  )

  expect_equal(
    angle_halfcircle_right()(theta + 2),
    c(-1 + pi, 0 + pi, 1, 2, 3, 4 - pi, 5 - pi)
  )

  expect_equal(
    angle_fixed()(theta),
    rep(0, length(theta))
  )

  expect_equal(
    angle_fixed(2)(theta),
    rep(2, length(theta))
  )

  expect_equal(
    angle_identity()(theta),
    theta
  )
})

test_that("none label placer", {
  lines <- list(
    "1" = list(x = c(.5, .75, .5, .25), y = c(.25, .5, .75, .5), id = rep(1, 4)),
    "2" = list(x = c(.5, 1, .5, 0), y = c(0, .5, 1, .5), id = rep(1, 4))
  )

  labels_data <- data.frame(
    index = 1:2,
    break_index = c(3, 7),
    break_id = c("1", "2"),
    label = c("a", "b"),
    stringsAsFactors = FALSE
  )

  lp <- label_placer_none()
  out <- lp(lines, labels_data)
  expect_equal(nrow(out), 0)
})


test_that("manual label placer", {
  lines <- list(
    "1" = list(x = c(.5, .75, .5, .25), y = c(.25, .5, .75, .5), id = rep(1, 4)),
    "2" = list(x = c(.5, 1, .5, 0), y = c(0, .5, 1, .5), id = rep(1, 4))
  )

  labels_data <- data.frame(
    index = 1:2,
    break_index = c(3, 7),
    break_id = c("1", "2"),
    label = c("a", "b"),
    stringsAsFactors = FALSE
  )

  lp <- label_placer_manual(
    breaks = c("1", "2", "3", "2"),
    x = 1:4,
    y = 4:1,
    theta = (1:4)-2
  )
  out <- lp(lines, labels_data)
  expect_equal(out$index, c(1, 2, 2))
  expect_equal(out$break_index, c(3, 7, 7))
  expect_equal(out$break_id, c("1", "2", "2"))
  expect_equal(out$label, c("a", "b", "b"))
  expect_equal(out$x, c(1, 2, 4))
  expect_equal(out$y, c(4, 3, 1))
  expect_equal(out$theta, c(-1, 0, 2))
})

# Two isolines id=1
test_that("middle label placer", {
  lines <- list(
    "1" = list(x = c(.5, .75, .5, .25), y = c(.25, .5, .75, .5), id = rep(1, 4)),
    "2" = list(x = c(.5, 1, .5, 0), y = c(0, .5, 1, .5), id = rep(1, 4))
  )

  labels_data <- data.frame(
    index = 1:2,
    break_index = c(3, 7),
    break_id = c("1", "2"),
    label = c("a", "b"),
    stringsAsFactors = FALSE
  )

  lp <- label_placer_middle()
  out <- lp(lines, labels_data)
  expect_equal(out$index, 1:2)
  expect_equal(out$break_index, c(3, 7))
  expect_equal(out$break_id, c("1", "2"))
  expect_equal(out$label, c("a", "b"))
  expect_equal(out$x, c(0.75, 1))
  expect_equal(out$y, c(0.5, 0.5))
  expect_equal(out$theta, c(pi / 2, pi / 2))
})

# One isoline with two id values (1 and 2)
test_that("middle label placer", {
  lines <- list(
    "1" = list(x = c(.5, .75, .5, .25, .5, 1, .5, 0), y = c(.25, .5, .75, .5, 0, .5, 1, .5), id = c(rep(1, 4), rep(2, 4)))
  )

  labels_data <- data.frame(
    index = 1,
    break_index = 3,
    break_id = "1",
    label = "a",
    stringsAsFactors = FALSE
  )

  lp <- label_placer_middle()
  out <- lp(lines, labels_data)
  expect_equal(out$index, c(1, 1))
  expect_equal(out$break_index, c(3, 3))
  expect_equal(out$break_id, c("1", "1"))
  expect_equal(out$label, c("a", "a"))
  expect_equal(out$x, c(0.75, 1))
  expect_equal(out$y, c(0.5, 0.5))
  expect_equal(out$theta, c(pi / 2, pi / 2))
})
