# Generic label placement function

The simple label placer processes separate isolines independently and
places labels for each line using a placer function that does the actual
placement work. This label placer is not meant to be used by end users,
but rather facilitates the development of new label placers, such as
[`label_placer_minmax()`](http://isoband.r-lib.org/reference/label_placer.md).

## Usage

``` r
label_placer_simple(lines, labels_data, placer_fun)
```

## Arguments

- lines:

  Isolines object for which labels should be placed.

- labels_data:

  A data frame containing information about which labels should be
  placed.

- placer_fun:

  A function that takes an individual isoline plus its associated break
  id as input and returns a data frame specifying label positions. The
  data frame should have three columns called `x`, `y`, and `theta`. `x`
  and `y` specify the label position, and `theta` specifies the label
  angle in radians. The data frame can have multiple rows, which results
  in the same label being placed in multiple locations.
