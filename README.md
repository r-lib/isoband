
<!-- README.md is generated from README.Rmd. Please edit that file -->

# isoband <img src="man/figures/isoband-logo.png" align="right" style="padding-left:10px;background-color:white;width:120px" />

[![CRAN
status](https://www.r-pkg.org/badges/version/isoband)](https://cran.r-project.org/package=isoband)
[![Build
Status](https://travis-ci.org/clauswilke/isoband.svg?branch=master)](https://travis-ci.org/clauswilke/isoband)
[![Coverage
Status](https://img.shields.io/codecov/c/github/clauswilke/isoband/master.svg)](https://codecov.io/github/clauswilke/isoband?branch=master)

Generate contour lines (isolines) and contour polygons (isobands) from
regularly spaced grids containing elevation data.

## Installation

Install the latest official release from CRAN via:

``` r
install.packages("isoband")
```

Install the current development from github via:

``` r
remotes::install_github("clauswilke/isoband")
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

The isolining and isobanding algorithms have no problem with larger
datasets. Let’s calculate isolines and isobands for a photographic
image, convert to sf, and plot with ggplot2.

``` r
library(magick)
#> Linking to ImageMagick 6.9.9.39
#> Enabled features: cairo, fontconfig, freetype, lcms, pango, rsvg, webp
#> Disabled features: fftw, ghostscript, x11
library(sf)
#> Linking to GEOS 3.6.1, GDAL 2.1.3, PROJ 4.9.3
library(ggplot2)

sf_from_image <- function(image) {
  image_gray <- image %>% image_quantize(colorspace = "gray")
  image_raster <- as.raster(image_gray)
  d <- dim(image_raster)
  m <- matrix(c((255-col2rgb(image_raster)[1,])), nrow = d[1], ncol = d[2], byrow = TRUE)
  b <- isobands(1:d[2], d[1]:1, m, 20*(0:13), 20*(1:14))
  bands <- iso_to_sfg(b)
  data <- st_sf(
    level = letters[1:length(bands)],
    geometry = st_sfc(bands)
  )
}

img <- image_resize(image_read(system.file("extdata", "ocean-cat.jpg", package = "isoband")), "200x200")
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
  p2 + scale_fill_viridis_d(
    aesthetics = c("color", "fill"), option = "B", guide = "none",
    direction = -1
  ),
  p2 + scale_fill_viridis_d(
    aesthetics = c("color", "fill"), option = "D", guide = "none",
    direction = -1
  ),
  p3,
  scale = 0.9
)
```

![](man/figures/README-polygon-cat-1.png)<!-- -->
