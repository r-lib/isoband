#' Efficient calculation of isolines and isobands from elevation grid
#'
#' @param x Numeric vector specifying the x locations of the grid points.
#' @param y Numeric vector specifying the y locations of the grid points.
#' @param z Numeric matrix specifying the elevation values for each grid point.
#' @param levels_low,levels_high Numeric vectors of minimum/maximum z values
#'   for which isobands should be generated.
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
  stats::setNames(out, paste0(levels_low, "-", levels_high))
}

#' @rdname isobands
#' @param levels Numeric vector of z values for which isolines should be generated.
#' @export
isolines <- function(x, y, z, levels) {
  out <- isolines_impl(x, y, z, levels)
  stats::setNames(out, levels)
}
