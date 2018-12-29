Rcpp::sourceCpp("~/Dropbox/projects/R_programming/contouring/contour_lines.cpp")

library(ggplot2)

contour_lines <- function(m, level = NULL) {
  if (is.null(level)) {
    level <- scales::pretty_breaks(10)(m)
  }
  purrr::map_dfr(level, ~single_contour_lines(1:ncol(m), nrow(m):1, m, .x))
}

# example of contour lines like current `geom_contour()`

df <- contour_lines(volcano)
ggplot(df, aes(x = x0, y = y0, xend = x1, yend = y1, color = level)) + 
  geom_segment() +
  scale_color_viridis_c()

# example of two contour bands, one from 120 to 140 and one from 150 to 152

df1 <- single_contour_bands(1:ncol(volcano), nrow(volcano):1, volcano, 120, 140)
df2 <- single_contour_bands(1:ncol(volcano), nrow(volcano):1, volcano, 150, 152)
df3 <- contour_lines(volcano, level = c(120, 140, 150, 152))
ggplot(mapping = aes(x, y, group = id)) + 
  geom_polygon(data = df1, color = "lightblue", fill = "lightblue") +
  geom_polygon(data = df2, color = "tomato", fill = "tomato") +
  geom_segment(data = df3, aes(x = x0, y = y0, xend = x1, yend = y1), size = 0.2, inherit.aes = FALSE)

# benchmark comparison to grDevices::contourLines()

microbenchmark::microbenchmark(
  grDevices::contourLines(1:ncol(volcano), 1:nrow(volcano), volcano, levels = 120),
  single_contour_lines(1:ncol(volcano), 1:nrow(volcano), volcano, 120),
  single_contour_bands(1:ncol(volcano), 1:nrow(volcano), volcano, 120, 140)
)
