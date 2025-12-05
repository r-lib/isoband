# Standardize label angles

Function factories that return functions to standardize rotation angles
to specific angle ranges.

## Usage

``` r
angle_halfcircle_bottom()

angle_halfcircle_right()

angle_fixed(theta = 0)

angle_identity()
```

## Arguments

- theta:

  Fixed angle, in radians.

## Details

`angle_halfcircle_bottom()` standardizes angles to (-pi/2, pi/2\].

`angle_halfcircle_right()` standardizes angles to (0, pi\].

`angle_fixed()` sets all angles to a fixed value (0 by default).

`angle_identity()` does not modify any angles.
