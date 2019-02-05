#' Render labeled isolines
#'
#' This function generates a grid grob that represents labeled isolines.
#'
#' @param lines Isolines, as produced by the [`isolines()`] function.
#' @param gp Grid graphical parameters.
#' @examples
#' x <- (1:ncol(volcano))/(ncol(volcano)+1)
#' y <- (nrow(volcano):1)/(nrow(volcano)+1)
#' lines <- isolines(x, y, volcano, 10*(10:18))
#'
#' library(grid)
#' grid.newpage()
#'
#' # make some colored background
#' grid.rect(width = .7, height = .7, gp = gpar(fill = "cornsilk1", col = NA))
#' grid.rect(width = .5, height = .5, gp = gpar(fill = "cornsilk2", col = NA))
#' grid.rect(width = .3, height = .3, gp = gpar(fill = "cornsilk3", col = NA))
#'
#' # draw labeled lines
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

  # calculate label width and height in npc units
  label_width <- convertWidth(stringWidth(label), "npc", valueOnly = TRUE)
  label_height <- convertHeight(stringHeight(label) + stringDescent(label), "npc", valueOnly = TRUE)

  # find minimum in contour line
  idx <- which(data$y == min(data$y))[1]

  # calculate label position
  pos <- label_position(data, idx, 1)
  clipped <- clip_lines(
    data$x, data$y, data$id, pos$center,
    label_width, label_height, pos$theta
  )

  rot <- 360*pos$theta/(2*pi)

  if (rot <= -90) {
    rot <- 180 + rot
  } else if (rot > 90) {
    rot <- rot - 180
  }

  grobTree(
    polylineGrob(clipped$x, clipped$y, clipped$id),
    textGrob(label, pos$center[1], pos$center[2], rot = rot)
  )
}

# Calculate the position and rotation of a label based
# on the x,y data and the index position in the data.
# The variable `n` sets the neighborhood size, n = 2 means
# two points in either direction are used.
label_position <- function(data, idx, n = 2) {
  imin <- max(idx - n, 1)
  imax <- min(idx + n, length(data$x))
  x <- data$x[imin:imax]
  y <- data$y[imin:imax]
  xave <- mean(x)
  yave <- mean(y)
  m <- cbind(x - xave, y - yave)
  v <- svd(m)$v
  list(center = c(xave, yave), theta = atan2(v[2], v[1]))
}
