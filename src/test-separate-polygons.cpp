#include <testthat.h>

#include "polygon.h"
#include "separate-polygons.h"

/* Things to regression-test for point_in_polygon():
 *  Degenerate polygon, single point
 *  Point aligned with point of one segment, inside
 *  Point aligned with point of one segment, outside
 */


context("Point in polygon") {
  test_that("Simple square") {
    polygon poly = {
      point(0, 0),
      point(0, 1),
      point(1, 1),
      point(1, 0),
      point(0, 0)
    };

    expect_true(point_in_polygon(point(0.5, 0.5), poly) == inside);
    expect_true(point_in_polygon(point(-0.5, 0.5), poly) == outside);
    expect_true(point_in_polygon(point(1.5, 0.5), poly) == outside);
    expect_true(point_in_polygon(point(0.5, -0.5), poly) == outside);
    expect_true(point_in_polygon(point(0.5, 1.5), poly) == outside);
    expect_true(point_in_polygon(point(-1, 1), poly) == outside);
    expect_true(point_in_polygon(point(2, 1), poly) == outside);
    expect_true(point_in_polygon(point(-1, 0), poly) == outside);
    expect_true(point_in_polygon(point(2, 0), poly) == outside);
    expect_true(point_in_polygon(point(0, 0), poly) == undetermined);
    expect_true(point_in_polygon(point(1, 0), poly) == undetermined);
    expect_true(point_in_polygon(point(0, 1), poly) == undetermined);
    expect_true(point_in_polygon(point(1, 1), poly) == undetermined);
  }

  test_that("Simple diamond") {
    polygon poly = {
      point(0, -.5),
      point(-.5, 0),
      point(0, .5),
      point(.5, 0),
      point(0, -.5)
    };

    expect_true(point_in_polygon(point(0, 0), poly) == inside);
    expect_true(point_in_polygon(point(-1, 0), poly) == outside);
    expect_true(point_in_polygon(point(1, 0), poly) == outside);
    expect_true(point_in_polygon(point(-.3, -.3), poly) == outside);
    expect_true(point_in_polygon(point(-.3, .3), poly) == outside);
    expect_true(point_in_polygon(point(.3, .3), poly) == outside);
    expect_true(point_in_polygon(point(.3, -.3), poly) == outside);
    expect_true(point_in_polygon(point(-.2, -.2), poly) == inside);
    expect_true(point_in_polygon(point(-.2, .2), poly) == inside);
    expect_true(point_in_polygon(point(.2, .2), poly) == inside);
    expect_true(point_in_polygon(point(.2, -.2), poly) == inside);
    expect_true(point_in_polygon(point(0, -.5), poly) == undetermined);
    expect_true(point_in_polygon(point(-.5, 0), poly) == undetermined);
    expect_true(point_in_polygon(point(0, .5), poly) == undetermined);
    expect_true(point_in_polygon(point(.5, 0), poly) == undetermined);
    expect_true(point_in_polygon(point(-.25, -.25), poly) == undetermined);
    expect_true(point_in_polygon(point(-.25, .25), poly) == undetermined);
    expect_true(point_in_polygon(point(.25, .25), poly) == undetermined);
    expect_true(point_in_polygon(point(.25, -.25), poly) == undetermined);
  }

  test_that("Degenerate polygon: horizontal line") {
    polygon poly = {
      point(0, 0),
      point(1, 0),
      point(2, 0),
      point(0, 0)
    };

    expect_true(point_in_polygon(point(-.5, 0), poly) == outside);
    expect_true(point_in_polygon(point(2.5, 0), poly) == outside);
    expect_true(point_in_polygon(point(0.5, 0), poly) == undetermined);
    expect_true(point_in_polygon(point(1.5, 0), poly) == undetermined);
  }

  test_that("Degenerate polygon: vertical line") {
    polygon poly = {
      point(.5, 2),
      point(.5, 1),
      point(.5, .5),
      point(.5, 2)
    };

    expect_true(point_in_polygon(point(0, 2), poly) == outside);
    expect_true(point_in_polygon(point(0, 1.5), poly) == outside);
    expect_true(point_in_polygon(point(0, 1), poly) == outside);
    expect_true(point_in_polygon(point(0, .8), poly) == outside);
    expect_true(point_in_polygon(point(0, .5), poly) == outside);
    expect_true(point_in_polygon(point(0, .4), poly) == outside);
    expect_true(point_in_polygon(point(1, 1), poly) == outside);
    expect_true(point_in_polygon(point(.5, 2), poly) == undetermined);
    expect_true(point_in_polygon(point(.5, 1.5), poly) == undetermined);
    expect_true(point_in_polygon(point(.5, 1), poly) == undetermined);
    expect_true(point_in_polygon(point(.5, .5), poly) == undetermined);
  }
}

