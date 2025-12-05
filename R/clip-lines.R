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
#' @keywords internal
#' @export
clip_lines <- function(x, y, id, clip_boxes, asp = 1) {
  out = list(x = x, y = y, id = id)
  for (i in 1:nrow(clip_boxes)) {
    box <- clip_boxes[i, ]
    out <- clip_lines_impl(
      as.double(out$x),
      as.double(out$y),
      as.integer(out$id),
      as.double(box$x),
      as.double(box$y),
      as.double(box$width),
      as.double(box$height),
      as.double(box$theta),
      as.double(asp)
    )
  }

  out
}
