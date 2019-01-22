
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

## Example

We begin with a few simple examples. The two main workhorses of the
package are the functions `isolines()` and `isobands()`, respectively.
They return a list of isolines/isobands for each isolevel specified.
Each isoline/isoband consists of vectors of x and y coordinates, as well
as a vector of ids specifying which sets of coordinates should be
connected. This format can be handed directly to
`grid.polyline()`/`grid.path()` for drawing (see below).

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

isobands(1:ncol(m), 1:nrow(m), m, 0.5, 1.5)
#> $`0.5-1.5`
#> $`0.5-1.5`$x
#>  [1] 2.50 2.00 1.50 1.50 1.50 2.00 3.00 4.00 4.50 4.00 3.75 4.00 4.50 4.00
#> [15] 3.50 3.00 3.00 3.25 3.50 3.00 2.50 2.50
#> 
#> $`0.5-1.5`$y
#>  [1] 4.00 4.50 4.00 3.00 2.00 1.50 1.25 1.50 2.00 2.50 3.00 3.50 4.00 4.50
#> [15] 4.00 3.75 3.25 3.00 2.00 1.75 2.00 3.00
#> 
#> $`0.5-1.5`$id
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2
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

The algorithm has no problem with larger datasets. Let’s calculate a few
isolines/isobands for the volcano dataset.

``` r
m <- volcano
b <- isobands((1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1), m, c(80, 120, 150), c(110, 140, 152))
l <- isolines((1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1), m, c(110, 120, 140, 150, 152))

library(grid)
grid.newpage()
grid.path(b[[1]]$x, b[[1]]$y, b[[1]]$id, gp = gpar(fill = "cornsilk", col = NA))
grid.path(b[[2]]$x, b[[2]]$y, b[[2]]$id, gp = gpar(fill = "lightblue", col = NA))
grid.path(b[[3]]$x, b[[3]]$y, b[[3]]$id, gp = gpar(fill = "tomato", col = NA))
for (i in seq_along(l)) {
  grid.polyline(l[[i]]$x, l[[i]]$y, l[[i]]$id)
}
```

<img src="man/figures/README-volcano-1.png" width="75%" />

Isolining is about as fast as `grDevices::contourLines()`, isobanding is
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
#>  1.629629 1.902334 2.916969 2.435455 3.156160 17.337077   100
#>  1.706788 1.859153 2.240307 2.185938 2.508821  3.658807   100
#>  4.356943 4.723402 5.443805 5.089021 6.081608  8.325564   100
```
