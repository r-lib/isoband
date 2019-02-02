#include <testthat.h>

#include "polygon.h"
#include "labels.h"

#include <cmath>

bool near_equal(double x1, double x2) {
  return (std::fabs(x1 - x2) < 1e-9);
}

context("Crop to unit box") {
  test_that("One point inside") {
    point crop1, crop2;
    segment_crop_type result;

    result = crop_to_unit_box(point(-.5, .5), point(.5, .5), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .5));
  }
}
