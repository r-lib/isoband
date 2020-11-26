test_that("basic functions", {
  m <- matrix(c(0, 0, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 1, 2, 1, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 0, 0), 5, 5, byrow = TRUE)

  l <- isolines((1:5)/6, (5:1)/6, m, c(.5, 1.5))

  # incorrect number of labels
  expect_error(
    isolines_grob(l, labels = c("a", "b", "c")),
    "Number of labels must match the number of breaks."
  )

  # incorrect margin specification
  expect_error(
    isolines_grob(l, margin = 1:4),
    "must be a unit object of length four"
  )
  expect_error(
    isolines_grob(l, margin = grid::unit(1:3, "pt")),
    "must be a unit object of length four"
  )

  # default settings create two labels per line for this dataset
  g <- isolines_grob(l, label_placer = label_placer_minmax(n = 0))
  expect_equal(g$labels_data$break_id, c("0.5", "0.5", "1.5", "1.5"))
  expect_equal(g$labels_data$label, c("0.5", "0.5", "1.5", "1.5"))
  expect_equal(g$labels_data$x, rep(0.5, 4))
  expect_equal(g$labels_data$theta, rep(0, 4))
  expect_true(all(abs(g$labels_data$y - c(0.75, 0.25, 0.5833333, 0.4166667)) < 1e-7))

  # extra breaks are ignored
  g <- isolines_grob(l, breaks = c("0.5", "1.5", "2.5"), labels = c("a", "b", "c"))
  expect_equal(g$labels_data$break_id, c("0.5", "0.5", "1.5", "1.5"))
  expect_equal(g$labels_data$label, c("a", "a", "b", "b"))
})
