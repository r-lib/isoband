#' Generate 2d kernel density estimates
#'
#' Generate 2d kernel density estimates
#' @inheritParams ggplot2::stat_identity
#' @inheritParams stat_isolevels
#' @param n Number of grid points in each direction.
#' @param h Bandwidth (vector of length two). If `NULL`, estimated
#'   using [MASS::bandwidth.nrd()].
#' @section Computed variables:
#' Same as [`stat_isolevels()`], with the addition of:
#' \describe{
#'   \item{`density`}{the density estimate}
#'   \item{`ndensity`}{density estimate scaled to maximum of 1}
#'   \item{`z`}{the density estimate, identical to `density`}
#' }
#' @examples
#' # default uses `geom_tile()` for drawing
#' ggplot(faithful, aes(eruptions, waiting)) +
#'   stat_density_grid(aes(fill = stat(density))) +
#'   geom_point(size = 0.3, color = "white")
#'
#' # discretized colors
#' ggplot(faithful, aes(eruptions, waiting)) +
#'   stat_density_grid(aes(fill = stat(zmin))) +
#'   geom_point(size = 0.3, color = "white")
#'
#' # points
#' ggplot(faithful, aes(eruptions, waiting)) +
#'   stat_density_grid(
#'     aes(color = stat(density), size = stat(density)),
#'     geom = "point", n = 20, stroke = 0
#'   ) +
#'   geom_point(size = 0.75, shape = 21, fill = "white", color = "black")
#'
#' # contour bands
#' ggplot(faithful, aes(eruptions, waiting)) +
#'   stat_density_grid(
#'     aes(fill = stat(zmin)), geom = "isobands",
#'     color = "gray70"
#'   ) +
#'   geom_point(size = 0.3, color = "white")
#' @export
stat_density_grid <- function(mapping = NULL, data = NULL,
                              geom = "tile", position = "identity",
                              ...,
                              n = 100, h = NULL,
                              bins = NULL, binwidth = NULL, breaks = NULL,
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatDensitygrid,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      n = n,
      h = h,
      bins = NULL, binwidth = NULL, breaks = NULL,
      ...
    )
  )
}

#' @rdname stat_density_grid
#' @usage NULL
#' @export
stat_densitygrid <- stat_density_grid


#' @rdname stat_density_grid
#' @format NULL
#' @usage NULL
#' @export
StatDensitygrid <- ggproto("StatDensitygrid", Stat,
  required_aes = c("x", "y"),

  compute_group = function(data, scales, na.rm = FALSE, h = NULL,
                           n = 100, bins = NULL, binwidth = NULL, breaks = NULL) {
    if (is.null(h)) {
      h <- c(MASS::bandwidth.nrd(data$x), MASS::bandwidth.nrd(data$y))
    }

    dens <- MASS::kde2d(
      data$x, data$y, h = h, n = n,
      lims = c(scales$x$dimension(), scales$y$dimension())
    )
    df <- expand.grid(x = dens$x, y = dens$y)
    df$density <- as.vector(dens$z)
    df$ndensity <- df$density / max(df$density)
    df$z <- df$density
    df$group <- data$group[1]

    StatIsolevels$compute_group(df, scales, bins, binwidth, breaks)
  }
)
