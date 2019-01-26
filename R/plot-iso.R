#' Visualize a single isoband
#'
#' This function visualizes a single isoband calculated from a matrix. It is mainly useful
#' for debugging and visualizing the isobanding algorithm. See [`isobands()`] for more
#' examples.
#'
#' @param m input matrix
#' @param vlo lower cutoff for isobanding
#' @param vhi higher cutoff for isobanding
#' @param col_lo line color for lower cutoff
#' @param col_hi line color for higher cutoff
#' @param fill_lo fill color for points below the lower cutoff
#' @param fill_mid fill color for points between the two cutoffs
#' @param fill_hi fill color for points above the higher cutoff
#' @param fill_band fill color for the isoband
#' @param newpage boolean, indicating whether `grid.newpage()` should
#'   be called or not
#' @examples
#' m <- matrix(c(0, 0, 0, 0, 0, 0,
#'               0, 2, 2, 2, 2, 0,
#'               0, 2, 0, 0, 2, 0,
#'               0, 2, 0, 0, 2, 0,
#'               0, 2, 2, 2, 2, 0,
#'               0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)
#'
#' plot_iso(m, 0.5, 1.5)
#' @export
plot_iso <- function(m, vlo, vhi, fill_lo = "gray95", fill_mid = "gray50", fill_hi = "black",
                     fill_band = "cornsilk", col_lo = "black", col_hi = "black", newpage = TRUE) {
  x <- 0.05 + 0.9*(0:(ncol(m)-1))/(ncol(m)-1)
  y <- 0.05 + 0.9*((nrow(m)-1):0)/(nrow(m)-1)
  df_bands <- isobands(x, y, m, vlo, vhi)[[1]]
  df_lines_lo <- isolines(x, y, m, vlo)[[1]]
  df_lines_hi <- isolines(x, y, m, vhi)[[1]]
  df_points <- expand.grid(y = y, x = x)
  pfill <- c(ifelse(m < vlo, fill_lo, ifelse(m < vhi, fill_mid, fill_hi)))
  pcol <- c(ifelse(m < vlo, "black", ifelse(m < vhi, fill_mid, "black")))
  if (isTRUE(newpage)) grid.newpage()
  if (length(df_bands$x) > 0)
    grid.path(df_bands$x, df_bands$y, df_bands$id, gp = gpar(fill = fill_band, col = NA))
  if (length(df_lines_lo$x) > 0)
    grid.polyline(df_lines_lo$x, df_lines_lo$y, df_lines_lo$id,
                  gp = gpar(col = col_lo))
  if (length(df_lines_hi$x) > 0)
    grid.polyline(df_lines_hi$x, df_lines_hi$y, df_lines_hi$id,
                  gp = gpar(col = col_hi))
  grid.points(df_points$x, df_points$y, default.units = "npc", pch = 21, size = unit(0.5, "char"),
              gp = gpar(fill = pfill, col = pcol))
}
