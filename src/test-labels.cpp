#include <testthat.h>

#include "polygon.h"
#include "labels.h"

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
}
