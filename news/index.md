# Changelog

## isoband 0.3.0

- General upkeep
- Rewrite implementation to use cpp11

## isoband 0.2.7

CRAN release: 2022-12-20

- Remove compile-time dependency on testthat.

- Changed maintainer after the original author (Claus Wilke) donated
  codebase to r-lib.

### isoband 0.2.6

- Update to the latest testthat headers for compatibility with LLVM
  clang 15.0.0.

- Correct label angle for current aspect ratio in
  [`isolines_grob()`](http://isoband.r-lib.org/reference/isolines_grob.md)
  ([\#28](https://github.com/r-lib/isoband/issues/28),
  [@eliocamp](https://github.com/eliocamp)).

### isoband 0.2.5

- Add a new label placer function
  [`label_placer_middle()`](http://isoband.r-lib.org/reference/label_placer.md)
  ([\#24](https://github.com/r-lib/isoband/issues/24),
  [@jamarav](https://github.com/jamarav)).

- The vendored testthat/catch code now uses a constant value for the
  stack size rather than relying on `SIGSTKSZ`. See:
  <https://github.com/r-lib/testthat/issues/1373>

### isoband 0.2.4

- Remove testthat compile-time dependency.

### isoband 0.2.3

- Fix build for testthat 3.0.

### isoband 0.2.2

- Remove Rcpp dependency
  ([\#11](https://github.com/r-lib/isoband/issues/11),
  [@thomasp85](https://github.com/thomasp85)).

### isoband 0.2.1

- Improved clipping algorithm for
  [`clip_lines()`](http://isoband.r-lib.org/reference/clip_lines.md),
  less likely to experience numerical instabilities.

### isoband 0.2.0

- Added
  [`isolines_grob()`](http://isoband.r-lib.org/reference/isolines_grob.md)
  for drawing labeled isolines via the grid graphics system. A companion
  function
  [`isobands_grob()`](http://isoband.r-lib.org/reference/isobands_grob.md)
  is provided for convenience.

- Numerous minor fixes and improvements.

### isoband 0.1.0

First public release.
