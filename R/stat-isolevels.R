#' Create discrete isolevels
#'
#' @usage NULL
#' @format NULL
#' @export
StatIsolevels <- ggplot2::ggproto("StatIsolevels", ggplot2::Stat,
  required_aes = c("z"),

  compute_group = function(data, scales, n = 10) {
    breaks <- pretty(data$z, n)
    if (min(breaks) > min(data$z)) {
      breaks <- c(min(data$z), breaks)
    }
    if (max(breaks) < max(data$z)) {
      breaks <- c(breaks, max(data$z))
    }

    level <- cut(data$z, breaks, labels = FALSE)
    zmin <- breaks[level]
    zmax <- breaks[level + 1]

    cbind(data, level = level, zmin = zmin, zmax = zmax)
  }
)
