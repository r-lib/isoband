#include <testthat.h>

#include "polygon.h"
#include "clip-lines.h"

#include <cmath>

bool near_equal(double x1, double x2) {
  return (std::fabs(x1 - x2) < 1e-9);
}

context("Crop to unit box") {
  test_that("Both points inside") {
    point crop1, crop2;
    segment_crop_type result;

    result = crop_to_unit_box(point(.2, .3), point(.7, .6), crop1, crop2);
    expect_true(result == complete);
    result = crop_to_unit_box(point(.7, .6), point(.2, .3), crop1, crop2);
    expect_true(result == complete);
  }

  test_that("Both points trivially outside") {
    point crop1, crop2;
    segment_crop_type result;

    // to the left
    result = crop_to_unit_box(point(-.2, .3), point(-.5, 1.6), crop1, crop2);
    expect_true(result == none);
    result = crop_to_unit_box(point(-.2, .3), point(-.5, .6), crop1, crop2);
    expect_true(result == none);
    result = crop_to_unit_box(point(-.5, 1.6), point(-.2, -.3), crop1, crop2);
    expect_true(result == none);

    // to the right
    result = crop_to_unit_box(point(1.2, .3), point(1.5, 1.6), crop1, crop2);
    expect_true(result == none);
    result = crop_to_unit_box(point(1.2, .3), point(1.5, .6), crop1, crop2);
    expect_true(result == none);
    result = crop_to_unit_box(point(1.5, 1.6), point(1.2, -.3), crop1, crop2);
    expect_true(result == none);

    // above
    result = crop_to_unit_box(point(.3, 1.2), point(1.6, 1.5), crop1, crop2);
    expect_true(result == none);
    result = crop_to_unit_box(point(.3, 1.2), point(.6, 1.5), crop1, crop2);
    expect_true(result == none);
    result = crop_to_unit_box(point(1.6, 1.5), point(-.3, 1.2), crop1, crop2);
    expect_true(result == none);

    // below
    result = crop_to_unit_box(point(.3, -.2), point(1.6, -.5), crop1, crop2);
    expect_true(result == none);
    result = crop_to_unit_box(point(.3, -.2), point(.6, -.5), crop1, crop2);
    expect_true(result == none);
    result = crop_to_unit_box(point(1.6, -.5), point(-.3, -.2), crop1, crop2);
    expect_true(result == none);
  }

  test_that("One point inside") {
    point crop1, crop2;
    segment_crop_type result;

    // horizontal lines
    result = crop_to_unit_box(point(-.2, .5), point(.5, .5), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .5));

    result = crop_to_unit_box(point(1.2, .5), point(.5, .5), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, .5));

    // vertical lines
    result = crop_to_unit_box(point(.5, .5), point(.5, -.2), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, .5));
    expect_true(near_equal(crop1.y, 0));

    result = crop_to_unit_box(point(.5, .5), point(.5, 1.2), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, .5));
    expect_true(near_equal(crop1.y, 1));

    // diagonal lines through corners
    result = crop_to_unit_box(point(.5, .5), point(1.2, 1.2), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, 1));

    result = crop_to_unit_box(point(.5, .5), point(-.2, 1.2), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, 1));

    result = crop_to_unit_box(point(-.2, -.2), point(.5, .5), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, 0));

    result = crop_to_unit_box(point(1.2, -.2), point(.5, .5), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, 0));

    // top
    result = crop_to_unit_box(point(.2, .8), point(.5, 1.4), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, .3));
    expect_true(near_equal(crop1.y, 1));

    // top right
    result = crop_to_unit_box(point(.8, .8), point(1.1, 1.4), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, .9));
    expect_true(near_equal(crop1.y, 1));

    result = crop_to_unit_box(point(.8, .8), point(1.4, 1.1), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, .9));

    // right
    result = crop_to_unit_box(point(1.4, .5), point(.8, .2), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, .3));

    // bottom right
    result = crop_to_unit_box(point(1.4, -.1), point(.8, .2), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, .1));

    result = crop_to_unit_box(point(1.4, -.35), point(.8, .05), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, .875));
    expect_true(near_equal(crop1.y, 0));

    // bottom
    result = crop_to_unit_box(point(.2, .2), point(.5, -.4), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, .3));
    expect_true(near_equal(crop1.y, 0));

    // bottom left
    result = crop_to_unit_box(point(.2, .2), point(-.4, -.1), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .1));

    result = crop_to_unit_box(point(.2, .05), point(-.4, -.35), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, .125));
    expect_true(near_equal(crop1.y, 0));

    // left
    result = crop_to_unit_box(point(-.4, .5), point(.2, .2), crop1, crop2);
    expect_true(result == at_end);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .3));

    // top left
    result = crop_to_unit_box(point(.2, .8), point(-.1, 1.4), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, .1));
    expect_true(near_equal(crop1.y, 1));

    result = crop_to_unit_box(point(.2, .8), point(-.4, 1.1), crop1, crop2);
    expect_true(result == at_beginning);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .9));
  }

  test_that("Double intersections") {
    point crop1, crop2;
    segment_crop_type result;

    // horizontal lines
    result = crop_to_unit_box(point(-1, .5), point(2, .5), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .5));
    expect_true(near_equal(crop2.x, 1));
    expect_true(near_equal(crop2.y, .5));

    result = crop_to_unit_box(point(2, .5), point(-1, .5), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, .5));
    expect_true(near_equal(crop2.x, 0));
    expect_true(near_equal(crop2.y, .5));

    // vertical lines
    result = crop_to_unit_box(point(.5, -1), point(.5, 2), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, .5));
    expect_true(near_equal(crop1.y, 0));
    expect_true(near_equal(crop2.x, .5));
    expect_true(near_equal(crop2.y, 1));

    result = crop_to_unit_box(point(.5, 2), point(.5, -1), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, .5));
    expect_true(near_equal(crop1.y, 1));
    expect_true(near_equal(crop2.x, .5));
    expect_true(near_equal(crop2.y, 0));

    // diagonals through corner points
    result = crop_to_unit_box(point(-3, -3), point(2, 2), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, 0));
    expect_true(near_equal(crop2.x, 1));
    expect_true(near_equal(crop2.y, 1));

    result = crop_to_unit_box(point(-1, 2), point(3, -2), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, 1));
    expect_true(near_equal(crop2.x, 1));
    expect_true(near_equal(crop2.y, 0));

    // top left corner
    result = crop_to_unit_box(point(-.4, .4), point(.4, 1.2), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .8));
    expect_true(near_equal(crop2.x, .2));
    expect_true(near_equal(crop2.y, 1));

    result = crop_to_unit_box(point(.4, 1.2), point(-.4, .4), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, .2));
    expect_true(near_equal(crop1.y, 1));
    expect_true(near_equal(crop2.x, 0));
    expect_true(near_equal(crop2.y, .8));

    // top right corner
    result = crop_to_unit_box(point(1.4, .4), point(.6, 1.2), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, .8));
    expect_true(near_equal(crop2.x, .8));
    expect_true(near_equal(crop2.y, 1));

    result = crop_to_unit_box(point(.6, 1.2), point(1.4, .4), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, .8));
    expect_true(near_equal(crop1.y, 1));
    expect_true(near_equal(crop2.x, 1));
    expect_true(near_equal(crop2.y, .8));

    // bottom left corner
    result = crop_to_unit_box(point(-.4, .6), point(.4, -.2), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .2));
    expect_true(near_equal(crop2.x, .2));
    expect_true(near_equal(crop2.y, 0));

    result = crop_to_unit_box(point(.4, -.2), point(-.4, .6), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, .2));
    expect_true(near_equal(crop1.y, 0));
    expect_true(near_equal(crop2.x, 0));
    expect_true(near_equal(crop2.y, .2));

    // bottom right corner
    result = crop_to_unit_box(point(.4, -.4), point(1.2, .4), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, .8));
    expect_true(near_equal(crop1.y, 0));
    expect_true(near_equal(crop2.x, 1));
    expect_true(near_equal(crop2.y, .2));

    result = crop_to_unit_box(point(1.2, .4), point(.4, -.4), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, .2));
    expect_true(near_equal(crop2.x, .8));
    expect_true(near_equal(crop2.y, 0));

    // horizontally across
    result = crop_to_unit_box(point(-1, -.2), point(3, 1.4), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 0));
    expect_true(near_equal(crop1.y, .2));
    expect_true(near_equal(crop2.x, 1));
    expect_true(near_equal(crop2.y, .6));

    result = crop_to_unit_box(point(3, 1.4), point(-1, -.2), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, 1));
    expect_true(near_equal(crop1.y, .6));
    expect_true(near_equal(crop2.x, 0));
    expect_true(near_equal(crop2.y, .2));

    // vertically across
    result = crop_to_unit_box(point(-.2, -1), point(1.4, 3), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, .2));
    expect_true(near_equal(crop1.y, 0));
    expect_true(near_equal(crop2.x, .6));
    expect_true(near_equal(crop2.y, 1));

    result = crop_to_unit_box(point(1.4, 3), point(-.2, -1), crop1, crop2);
    expect_true(result == in_middle);
    expect_true(near_equal(crop1.x, .6));
    expect_true(near_equal(crop1.y, 1));
    expect_true(near_equal(crop2.x, .2));
    expect_true(near_equal(crop2.y, 0));
  }

  test_that("Points non-trivially outside") {
    point crop1, crop2;
    segment_crop_type result;

    result = crop_to_unit_box(point(-.2, .9), point(.1, 1.2), crop1, crop2);
    expect_true(result == none);

    result = crop_to_unit_box(point(1.2, .9), point(.9, 1.2), crop1, crop2);
    expect_true(result == none);

    result = crop_to_unit_box(point(-.2, .1), point(.1, -.2), crop1, crop2);
    expect_true(result == none);

    result = crop_to_unit_box(point(1.2, .1), point(.9, -.2), crop1, crop2);
    expect_true(result == none);
  }
}


context("Transform to unit box") {
  test_that("Simple transformations work") {
    unitbox_transformer t(point(1, 1), point(2, 2), point(0, 2));

    point p = t.transform(point(1, 2));
    point p2 = t.inv_transform(p);
    expect_true(near_equal(p.x, .5));
    expect_true(near_equal(p.y, .5));
    expect_true(near_equal(p2.x, 1));
    expect_true(near_equal(p2.y, 2));

    p = t.transform(point(1, 3));
    p2 = t.inv_transform(p);
    expect_true(near_equal(p.x, 1));
    expect_true(near_equal(p.y, 1));
    expect_true(near_equal(p2.x, 1));
    expect_true(near_equal(p2.y, 3));
  }

  test_that("Transformations from/to rhomboid work") {
    unitbox_transformer t(point(1, 1), point(2, 1), point(2, 2));

    point p = t.transform(point(2, 2));
    point p2 = t.inv_transform(p);
    expect_true(near_equal(p.x, 0));
    expect_true(near_equal(p.y, 1));
    expect_true(near_equal(p2.x, 2));
    expect_true(near_equal(p2.y, 2));

    p = t.transform(point(3, 2));
    p2 = t.inv_transform(p);
    expect_true(near_equal(p.x, 1));
    expect_true(near_equal(p.y, 1));
    expect_true(near_equal(p2.x, 3));
    expect_true(near_equal(p2.y, 2));
  }

  /*
  // the following tests don't work properly, because `expect_error()`
  // doesn't catch calls to Rf_error(), it only catches exceptions.
  test_that("Singular transformations are caught") {
    expect_error(
      // box without width
      unitbox_transformer(point(1, 1), point(1, 1), point(0, 2))
    );

    expect_error(
      // box without height
      unitbox_transformer(point(1, 1), point(2, 2), point(1, 1))
    );

    expect_error(
      // singular inverse transform
      unitbox_transformer(point(1, 1), point(2, 2), point(2, 2))
    );
  }
*/
}
