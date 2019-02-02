#' Render labeled isolines
#'
#' This function generates a grid grob that represents labeled isolines.
#'
#' @param lines Isolines, as produced by the [`isolines()`] function.
#' @param gp Grid graphical parameters.
#' @examples
#' m <- matrix(c(0, 0, 0, 0, 0, 0,
#'               0, 0, 0, 1, 1, 0,
#'               0, 0, 1, 2, 1, 0,
#'               0, 1, 2, 0, 0, 0,
#'               0, 0, 0, 3, 0, 0,
#'               0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)
#'
#' lines <- isolines(
#'   (1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1),
#'   m, c(0.5, 1.5, 2.5, 3.5)
#' )
#'
#' library(grid)
#' grid.newpage()
#' grid.draw(isolines_grob(lines))
#' @export
isolines_grob <- function(lines, gp = gpar()) {
  gTree(lines = lines, gp = gp, cl = "isolines_grob")
}

makeContent.isolines_grob <- function(x) {
  grobs <- mapply(
    labeled_polyline_grob,
    x$lines,
    names(x$lines),
    SIMPLIFY = FALSE
  )

  children <- do.call(gList, grobs)
  setChildren(x, children)
}

labeled_polyline_grob <- function(data, label) {
  if (length(data$x) == 0) {
    return(NULL)
  }

  idx <- which(data$y == min(data$y))[1]

  grobTree(
    polylineGrob(data$x, data$y, data$id),
    textGrob(label, data$x[idx], data$y[idx])
  )
}



# Test code ---------------------------------------------------------------
# Some test code, to be deleted eventually
# Demonstrates how we can obtain the size of a text label relative to
# the size of the viewport in which it is drawn

test_grob <- function(label, gp = gpar()) {
  gTree(label = label, gp = gp, cl = "test_grob")
}

makeContent.test_grob <- function(x) {
  grob_width_in <- convertWidth(unit(1, "npc"), "in", valueOnly = TRUE)
  grob_height_in <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)

  label_width_in <- convertWidth(stringWidth(x$label), "in", valueOnly = TRUE)
  label_height_in <- convertHeight(stringHeight(x$label) + stringDescent(x$label), "in", valueOnly = TRUE)

  #print(glue::glue("grob: ({grob_width_in}, {grob_height_in}); label: ({label_width_in}, {label_height_in})"))
  textGrob(x$label)
}

#grid.newpage()
#grid.draw(test_grob("hello!"))
#grid.draw(test_grob("world", gp = gpar(fontsize = 24)))


