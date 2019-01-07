#' Bin a 2d elevation raster into bands of similar height
#'
#' Bin a 2d elevation raster into bands of similar height.
#' @inheritParams ggplot2::stat_identity
#' @param bins Number of bins for discretization. Has priority
#'   over `binwidth`.
#' @param binwidth Binwidth used during discretization.
#' @param breaks Explicit bin boundaries to use for discretization.
#'   Has priority over both `bins` and `binwidth`.
#' @section Computed variables:
#' \describe{
#'  \item{`level`}{height of contour band, expressed as integer}
#'  \item{`zmin`}{minimum z value of contour band}
#'  \item{`zmax`}{maximum z value of contour band}
#'  \item{`nzmax`}{maximum z value of contour band, normalized to a maximum of 1}
#' }
#' @export
stat_isolevels <- function(mapping = NULL, data = NULL,
                         geom = "contour", position = "identity",
                         ...,
                         bins = NULL, binwidth = NULL, breaks = NULL,
                         na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatContour,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      bins = bins,
      binwidth = binwidth,
      breaks = breaks,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname stat_isolevels
#' @usage NULL
#' @format NULL
#' @export
StatIsolevels <- ggplot2::ggproto("StatIsolevels", ggplot2::Stat,
  required_aes = c("z"),

  compute_group = function(data, scales, bins = NULL, binwidth = NULL,
                           breaks = NULL) {
    # expand range by 1% in each direction
    r <- range(data$z)
    d <- diff(r)
    if (d == 0) d <- 1
    er <- c(r[1] - .01*d, r[2] + .01*d)

    # if no parameters set, use pretty bins
    if (is.null(bins) && is.null(binwidth) && is.null(breaks)) {
      breaks <- pretty(er, 10)
    }
    # if provided, use bins to calculate binwidth
    if (!is.null(bins)) {
      binwidth <- diff(er) / bins
    }
    # if necessary, compute breaks from binwidth
    if (is.null(breaks)) {
      breaks <- scales::fullseq(er, binwidth)
    }

    level <- cut(data$z, breaks, labels = FALSE)
    zmin <- breaks[level]
    zmax <- breaks[level + 1]
    nzmax <- zmax/max(zmax)

    cbind(data, level = level, zmin = zmin, zmax = zmax, nzmax = nzmax)
  }
)
