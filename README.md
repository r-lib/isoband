
<!-- README.md is generated from README.Rmd. Please edit that file -->

# isoband

[![Build
Status](https://travis-ci.org/clauswilke/isoband.svg?branch=master)](https://travis-ci.org/clauswilke/isoband)
[![Coverage
Status](https://img.shields.io/codecov/c/github/clauswilke/isoband/master.svg)](https://codecov.io/github/clauswilke/isoband?branch=master)

Generate contour lines (isolines) and contour polygons (isobands) from
regularly spaced grids containing elevation data.

## Installation

Install from github with:

``` r
devtools::install_github("clauswilke/isoband")
```

## Examples

The two main workhorses of the package are the functions `isolines()`
and `isobands()`, respectively. They return a list of isolines/isobands
for each isolevel specified. Each isoline/isoband consists of vectors of
x and y coordinates, as well as a vector of ids specifying which sets of
coordinates should be connected. This format can be handed directly to
`grid.polyline()`/`grid.path()` for drawing. However, we can also
convert the output to spatial features and draw with ggplot2 (see
below).

``` r
library(isoband)

m <- matrix(c(0, 0, 0, 0, 0,
              0, 1, 2, 1, 0,
              0, 1, 2, 0, 0,
              0, 1, 0, 1, 0,
              0, 0, 0, 0, 0), 5, 5, byrow = TRUE)

isolines(1:ncol(m), 1:nrow(m), m, 0.5)
#> $`0.5`
#> $`0.5`$x
#>  [1] 4.00 3.50 3.00 2.50 2.00 1.50 1.50 1.50 2.00 3.00 4.00 4.50 4.00 3.75
#> [15] 4.00 4.50 4.00
#> 
#> $`0.5`$y
#>  [1] 4.50 4.00 3.75 4.00 4.50 4.00 3.00 2.00 1.50 1.25 1.50 2.00 2.50 3.00
#> [15] 3.50 4.00 4.50
#> 
#> $`0.5`$id
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> 
#> attr(,"class")
#> [1] "isolines" "iso"

isobands(1:ncol(m), 1:nrow(m), m, 0.5, 1.5)
#> $`0.5:1.5`
#> $`0.5:1.5`$x
#>  [1] 2.50 2.00 1.50 1.50 1.50 2.00 3.00 4.00 4.50 4.00 3.75 4.00 4.50 4.00
#> [15] 3.50 3.00 3.00 3.25 3.50 3.00 2.50 2.50
#> 
#> $`0.5:1.5`$y
#>  [1] 4.00 4.50 4.00 3.00 2.00 1.50 1.25 1.50 2.00 2.50 3.00 3.50 4.00 4.50
#> [15] 4.00 3.75 3.25 3.00 2.00 1.75 2.00 3.00
#> 
#> $`0.5:1.5`$id
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2
#> 
#> 
#> attr(,"class")
#> [1] "isobands" "iso"
```

The function `plot_iso()` is a convenience function for debugging and
testing.

``` r
plot_iso(m, 0.5, 1.5)
```

<img src="man/figures/README-basic-example-plot-1.png" width="50%" />

A few more simple examples. Missing values and disconnected areas are
all handled correctly.

``` r
m <- matrix(c(0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 1, 0,
              0, 0, 1, 1, 1, 0,
              0, 1, 1, 0, 0, 0,
              0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)
plot_iso(m, 0.5, 1.5)
```

<img src="man/figures/README-basic-plotting-1.png" width="50%" />

``` r

m <- matrix(c(NA, NA, NA, 0, 0, 0,
              NA, NA, NA, 1, 1, 0,
               0,  0,  1, 1, 1, 0,
               0,  1,  1, 0, 0, 0,
               0,  0,  0, 1, 0, 0,
               0,  0,  0, 0, 0, 0), 6, 6, byrow = TRUE)
plot_iso(m, 0.5, 1.5)
```

<img src="man/figures/README-basic-plotting-2.png" width="50%" />

``` r

m <- matrix(c(0, 0, 1, 1,
              0, 1, 1, 1,
              1, 1, 0, 0,
              0, 0, 0.8, 0), 4, 4, byrow = TRUE)
plot_iso(m, 0.5, 1.5)
```

<img src="man/figures/README-basic-plotting-3.png" width="50%" />

The algorithm has no problem with larger datasets. Let’s calculate
isolines and isobands for the volcano dataset, convert to sf, and plot
with ggplot2.

``` r
library(ggplot2)
library(sf)
#> Linking to GEOS 3.6.1, GDAL 2.1.3, PROJ 4.9.3

m <- volcano
b <- isobands((1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1), m, 10*(9:19), 10*(10:20))
l <- isolines((1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1), m, 10*(10:19))

bands <- iso_to_sfg(b)
data_bands <- st_sf(
  level = 1:length(bands),
  geometry = st_sfc(bands)
)
lines <- iso_to_sfg(l)
data_lines <- st_sf(
  level = 2:(length(lines)+1),
  geometry = st_sfc(lines)
)

ggplot() +
  geom_sf(data = data_bands, aes(fill = level), color = NA, alpha = 0.7) +
  geom_sf(data = data_lines, color = "black") +
  scale_fill_viridis_c(guide = "none") +
  coord_sf(expand = FALSE)
```

<img src="man/figures/README-volcano-1.png" width="75%" />

We can also use this approach to convert bitmap images into polygons and
plot with ggplot2.

``` r
library(magick)
#> Linking to ImageMagick 6.9.9.39
#> Enabled features: cairo, fontconfig, freetype, lcms, pango, rsvg, webp
#> Disabled features: fftw, ghostscript, x11
library(sf)
library(ggplot2)
library(isoband)
library(colorspace)

sf_from_image <- function(image) {
  image_gray <- image %>% image_quantize(colorspace = "gray")
  image_raster <- as.raster(image_gray)
  d <- dim(image_raster)
  m <- matrix(c((255-col2rgb(image_raster)[1,])), nrow = d[1], ncol = d[2], byrow = TRUE)
  b <- isobands(1:d[2], d[1]:1, m, 20*(-1:13), 20*(0:14))
  bands <- iso_to_sfg(b)
  data <- st_sf(
    level = letters[1:length(bands)],
    geometry = st_sfc(bands)
  )
}

img <- image_resize(image_read("https://pbs.twimg.com/profile_images/905186381995147264/7zKAG5sY_400x400.jpg"), "200x200")
img_sf <- sf_from_image(img)

p1 <- ggplot(img_sf) + 
  geom_sf(color = "gray10", fill = NA, size = 0.05) + 
  coord_sf(expand = FALSE) +
  theme_gray() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks.length = grid::unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0)
  )

p2 <- ggplot(img_sf) + 
  geom_sf(aes(fill = level, color = level)) + 
  coord_sf(expand = FALSE) +
  theme_void()

p3 <- ggplot(img_sf) + 
  geom_sf(aes(fill = level), color = "gray30", size = 0.1) + 
  scale_fill_hue(aesthetics = c("color", "fill"), guide = "none", direction = -1) +
  coord_sf(expand = FALSE) +
  theme_void()

cowplot::plot_grid(
  p1,
  p2 + scale_fill_discrete_sequential(aesthetics = c("color", "fill"), palette = "Inferno", guide = "none", rev = TRUE),
  p2 + scale_fill_discrete_sequential(aesthetics = c("color", "fill"), palette = "Viridis", guide = "none", rev = TRUE),
  p3,
  scale = 0.9
)
```

![](man/figures/README-polygon-hadley-1.png)<!-- -->

## Performance

The code is written in C++ and performance is generally good. Isolining
is about as fast as `grDevices::contourLines()`, isobanding is
approximately 2.5 times slower.

``` r
microbenchmark::microbenchmark(
  grDevices::contourLines(1:ncol(volcano), 1:nrow(volcano), volcano, levels = 10*(10:18)),
  isolines(1:ncol(volcano), 1:nrow(volcano), volcano, 10*(10:18)),
  isobands(1:ncol(volcano), 1:nrow(volcano), volcano, 10*(9:17), 10*(10:18))
)
#> Unit: milliseconds
#>                                                                                            expr
#>  grDevices::contourLines(1:ncol(volcano), 1:nrow(volcano), volcano,      levels = 10 * (10:18))
#>                               isolines(1:ncol(volcano), 1:nrow(volcano), volcano, 10 * (10:18))
#>             isobands(1:ncol(volcano), 1:nrow(volcano), volcano, 10 * (9:17),      10 * (10:18))
#>       min       lq     mean   median       uq       max neval
#>  1.675161 1.812835 2.694624 2.129540 2.744745 16.622916   100
#>  1.761252 1.850351 2.429487 1.959520 2.346297  9.678578   100
#>  4.376592 4.630960 5.593479 4.959469 5.679508 13.972505   100
```
