context("test-isolines-grob")

test_that("basic functions", {
  m <- matrix(c(0, 0, 0, 0, 0, 0,
                0, 1, 1, 1, 1, 0,
                0, 1, 2, 2, 1, 0,
                0, 1, 2, 2, 1, 0,
                0, 1, 1, 1, 1, 0,
                0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)

  l <- isolines((1:6)/7, (6:1)/7, m, c(.5, 1.5))

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
  g <- isolines_grob(l)
  expect_equal(g$labels_data$break_id, c("0.5", "0.5", "1.5", "1.5"))
  expect_equal(g$labels_data$label, c("0.5", "0.5", "1.5", "1.5"))
})
