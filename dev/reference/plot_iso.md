# Visualize a single isoband

This function visualizes a single isoband calculated from a matrix. It
is mainly useful for debugging and visualizing the isobanding algorithm.
See [`isobands()`](http://isoband.r-lib.org/dev/reference/isobands.md)
for more examples.

## Usage

``` r
plot_iso(
  m,
  vlo,
  vhi,
  fill_lo = "gray95",
  fill_mid = "gray50",
  fill_hi = "black",
  fill_band = "cornsilk",
  col_lo = "black",
  col_hi = "black",
  newpage = TRUE
)
```

## Arguments

- m:

  input matrix

- vlo:

  lower cutoff for isobanding

- vhi:

  higher cutoff for isobanding

- fill_lo:

  fill color for points below the lower cutoff

- fill_mid:

  fill color for points between the two cutoffs

- fill_hi:

  fill color for points above the higher cutoff

- fill_band:

  fill color for the isoband

- col_lo:

  line color for lower cutoff

- col_hi:

  line color for higher cutoff

- newpage:

  boolean, indicating whether
  [`grid.newpage()`](https://rdrr.io/r/grid/grid.newpage.html) should be
  called or not

## Examples

``` r
m <- matrix(c(0, 0, 0, 0, 0, 0,
              0, 2, 2, 2, 2, 0,
              0, 2, 0, 0, 2, 0,
              0, 2, 0, 0, 2, 0,
              0, 2, 2, 2, 2, 0,
              0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)

plot_iso(m, 0.5, 1.5)
```
