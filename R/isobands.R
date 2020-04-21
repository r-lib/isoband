#' Efficient calculation of isolines and isobands from elevation grid
#'
#' @param x Numeric vector specifying the x locations of the grid points.
#' @param y Numeric vector specifying the y locations of the grid points.
#' @param z Numeric matrix specifying the elevation values for each grid point.
#' @param levels_low,levels_high Numeric vectors of minimum/maximum z values
#'   for which isobands should be generated. Any z values that are exactly
#'   equal to a value in `levels_low` are considered part of the corresponding
#'   isoband, but any z values that are exactly equal to a value in `levels_high`
#'   are not considered part of the corresponding isoband. In other words, the
#'   intervals specifying isobands are closed at their lower boundary and open
#'   at their upper boundary.
#' @seealso
#' [`plot_iso`]
#' @examples
#' library(grid)
#'
#' #' # one simple connected shape
#' m <- matrix(c(0, 0, 0, 0, 0, 0,
#'               0, 0, 0, 1, 1, 0,
#'               0, 0, 1, 1, 1, 0,
#'               0, 1, 1, 0, 0, 0,
#'               0, 0, 0, 1, 0, 0,
#'               0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)
#'
#' df_bands <- isobands((1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1), m, 0.5, 1.5)[[1]]
#' df_lines <- isolines((1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1), m, 0.5)[[1]]
#' g <- expand.grid(x = (1:ncol(m))/(ncol(m)+1), y = (nrow(m):1)/(nrow(m)+1))
#' grid.newpage()
#' grid.points(g$x, g$y, default.units = "npc", pch = 19, size = unit(0.5, "char"))
#' grid.path(df_bands$x, df_bands$y, df_bands$id, gp = gpar(fill = "cornsilk", col = NA))
#' grid.polyline(df_lines$x, df_lines$y, df_lines$id)
#'
#' # a similar plot can be generated with the plot_iso() function,
#' # which is useful for exploring how the algorithm works
#' plot_iso(m, 0.5, 1.5)
#'
#' # NAs are ignored
#' m <- matrix(c(NA, NA, NA, 0, 0, 0,
#'               NA, NA, NA, 1, 1, 0,
#'                0,  0,  1, 1, 1, 0,
#'                0,  1,  1, 0, 0, 0,
#'                0,  0,  0, 1, 0, 0,
#'                0,  0,  0, 0, 0, 0), 6, 6, byrow = TRUE)
#' plot_iso(m, 0.5, 1.5)
#'
#' # two separate shapes
#' m <- matrix(c(0, 0, 1, 1,
#'               0, 1, 1, 1,
#'               1, 1, 0, 0,
#'               0, 0, 0.8, 0), 4, 4, byrow = TRUE)
#' plot_iso(m, 0.5, 1.5)
#'
#' # shape with hole
#' m <- matrix(c(0, 0, 0, 0, 0, 0,
#'               0, 1, 1, 1, 1, 0,
#'               0, 1, 2, 2, 1, 0,
#'               0, 1, 2, 2, 1, 0,
#'               0, 1, 1, 1, 1, 0,
#'               0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)
#' plot_iso(m, 0.5, 1.5)
#' @export
isobands <- function(x, y, z, levels_low, levels_high) {
  nlow <- length(levels_low)
  nhigh <- length(levels_high)
  nmax <- max(nlow, nhigh)

  if ((nlow != nmax && nlow != 1) || (nhigh != nmax && nhigh != 1)) {
    stop("Vectors specifying isoband levels must be of equal length or of length 1", call. = FALSE)
  }
  levels_low <- rep_len(levels_low, nmax)
  levels_high <- rep_len(levels_high, nmax)

  # swap high and low levels when they're given in the wrong order
  idx <- levels_high < levels_low
  if (any(idx)) {
    levels_tmp <- levels_high
    levels_high[idx] <- levels_low[idx]
    levels_low[idx] <- levels_tmp[idx]
  }

  out <- isobands_impl(x, y, z, levels_low, levels_high)
  structure(
    stats::setNames(out, paste0(levels_low, ":", levels_high)),
    class = c("isobands", "iso")
  )
}

#' @rdname isobands
#' @param levels Numeric vector of z values for which isolines should be generated.
#' @export
isolines <- function(x, y, z, levels) {
  out <- isolines_impl(x, y, z, levels)
  structure(
    stats::setNames(out, levels),
    class = c("isolines", "iso")
  )
}

isobands_impl <- function(x, y, z, value_low, value_high) {
  mode(z) <- "numeric"
  .Call(
    `isobands_impl_c`, as.numeric(x), as.numeric(y), z,
    as.numeric(value_low), as.numeric(value_high), PACKAGE = "isoband"
  )
}

isolines_impl <- function(x, y, z, value) {
  mode(z) <- "numeric"
  .Call(
    `isolines_impl_c`, as.numeric(x), as.numeric(y), z,
    as.numeric(value), PACKAGE = "isoband"
  )
}
