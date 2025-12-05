# Set up a label placement strategy

These functions set up various label placement strategies.

## Usage

``` r
label_placer_minmax(
  placement = "tb",
  rot_adjuster = angle_halfcircle_bottom(),
  n = 2
)

label_placer_none()

label_placer_manual(breaks, x, y, theta)

label_placer_middle(rot_adjuster = angle_halfcircle_bottom())
```

## Arguments

- placement:

  String consisting of any combination of the letters "t", "r", "b", "l"
  indicating the placement of labels at the top, to the right, at the
  bottom, to the left of the isoline.

- rot_adjuster:

  Function that standardizes the rotation angles of the labels. See e.g.
  [`angle_halfcircle_bottom()`](angle_halfcircle_bottom.md).

- n:

  Size of the point neighborhood over which the rotation angle should be
  calculated.

- breaks:

  Character vector specifying the isolines to be labeled, as in
  [`isolines_grob()`](isolines_grob.md).

- x, y, theta:

  Numeric vectors specifying the x and y positions and angles (in
  radians) for each label corresponding to each break.

## Details

`label_placer_minmax()` places labels at the horizontal or vertical
minima or maxima of the respective isolines.

`label_placer_none()` places no labels at all.

`label_placer_manual()` places labels at manually defined locations.

`label_placer_middle()` places labels at the middle of each isoline.
