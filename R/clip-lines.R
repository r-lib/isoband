#' Clip lines so they don't run into a set of boxes.
#'
#' Clip lines so they don't run into a set of boxes. Useful for labeling isolines,
#' as it allows removal of line segments that would run into any text labels.
#' @param x Numeric vector of x coordinates
#' @param y Numeric vector of y coordinates
#' @param id Integer vector of id numbers indicating which lines are connected
#' @param clip_boxes Data frame specifying the locations of boxes to clip to.
#'   Should have five columns, named `x`, `y`, `width`, `height`, `theta`, which
#'   specify the x and y positions of each box midpoint, as well as the box width,
#'   box height, and box angle in radians. Each box is specified by one data
#'   row.
#' @param asp Aspect ratio (width/height) of the target canvas. This is used to convert
#'   widths to heights and vice versa for rotated boxes
#' @export
clip_lines <- function(x, y, id, clip_boxes, asp = 1) {
  out = list(x = x, y = y, id = id)
  for (i in 1:nrow(clip_boxes)) {
    box <- clip_boxes[i, ]
    out <- clip_lines_impl(
      out$x, out$y, out$id,
      box$x, box$y, box$width, box$height, box$theta,
      asp
    )
  }

  out
}

clip_lines_impl <- function(x, y, id, p_mid_x, p_mid_y, width, height, theta, asp = 1) {
  .Call(
    `clip_lines_impl_c`, as.numeric(x), as.numeric(y), as.integer(id),
    as.numeric(p_mid_x), as.numeric(p_mid_y), as.numeric(width),
    as.numeric(height), as.numeric(theta), as.numeric(asp), PACKAGE = "isoband"
  )
}
