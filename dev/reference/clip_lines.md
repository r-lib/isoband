# Clip lines so they don't run into a set of boxes.

Clip lines so they don't run into a set of boxes. Useful for labeling
isolines, as it allows removal of line segments that would run into any
text labels.

## Usage

``` r
clip_lines(x, y, id, clip_boxes, asp = 1)
```

## Arguments

- x:

  Numeric vector of x coordinates

- y:

  Numeric vector of y coordinates

- id:

  Integer vector of id numbers indicating which lines are connected

- clip_boxes:

  Data frame specifying the locations of boxes to clip to. Should have
  five columns, named `x`, `y`, `width`, `height`, `theta`, which
  specify the x and y positions of each box midpoint, as well as the box
  width, box height, and box angle in radians. Each box is specified by
  one data row.

- asp:

  Aspect ratio (width/height) of the target canvas. This is used to
  convert widths to heights and vice versa for rotated boxes
