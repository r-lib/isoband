#include <testthat.h>

#include "polygon.h"
#include "separate-polygons.h"

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

    // alternative version
    polygon poly2 = {
      point(1, 0),
      point(2, 0),
      point(0, 0),
      point(1, 0)
    };

    expect_true(point_in_polygon(point(-.5, 0), poly2) == outside);
    expect_true(point_in_polygon(point(2.5, 0), poly2) == outside);
    expect_true(point_in_polygon(point(0.5, 0), poly2) == undetermined);
    expect_true(point_in_polygon(point(1.5, 0), poly2) == undetermined);
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

    // alternative version
    polygon poly2 = {
      point(.5, 1),
      point(.5, .5),
      point(.5, 2),
      point(.5, 1)
    };

    expect_true(point_in_polygon(point(0, 2), poly2) == outside);
    expect_true(point_in_polygon(point(0, 1.5), poly2) == outside);
    expect_true(point_in_polygon(point(0, 1), poly2) == outside);
    expect_true(point_in_polygon(point(0, .8), poly2) == outside);
    expect_true(point_in_polygon(point(0, .5), poly2) == outside);
    expect_true(point_in_polygon(point(0, .4), poly2) == outside);
    expect_true(point_in_polygon(point(1, 1), poly2) == outside);
    expect_true(point_in_polygon(point(.5, 2), poly2) == undetermined);
    expect_true(point_in_polygon(point(.5, 1.5), poly2) == undetermined);
    expect_true(point_in_polygon(point(.5, 1), poly2) == undetermined);
    expect_true(point_in_polygon(point(.5, .5), poly2) == undetermined);
  }

  test_that("Degenerate polygon: point") {
    polygon poly = {
      point(0, 0),
      point(0, 0)
    };

    expect_true(point_in_polygon(point(-1, 0), poly) == outside);
    expect_true(point_in_polygon(point(1, 0), poly) == outside);
    expect_true(point_in_polygon(point(0, -1), poly) == outside);
    expect_true(point_in_polygon(point(0, 1), poly) == outside);
    expect_true(point_in_polygon(point(0.5, 0.5), poly) == outside);
    expect_true(point_in_polygon(point(0, 0), poly) == undetermined);
  }


  test_that("Multiple flat line segments") {
    polygon poly = {
      point(0, 2),
      point(1, 1),
      point(2, 1),
      point(3, 1),
      point(4, 1),
      point(4, 0),
      point(0, 0),
      point(0, 2)
    };

    expect_true(point_in_polygon(point(-1, 1), poly) == outside);
    expect_true(point_in_polygon(point(5, 1), poly) == outside);
    expect_true(point_in_polygon(point(0.5, 1), poly) == inside);

    // alternative version
    polygon poly2 = {
      point(1, 1),
      point(2, 1),
      point(3, 1),
      point(4, 1),
      point(4, 0),
      point(0, 0),
      point(0, 2),
      point(1, 1)
    };

    expect_true(point_in_polygon(point(-1, 1), poly2) == outside);
    expect_true(point_in_polygon(point(5, 1), poly2) == outside);
    expect_true(point_in_polygon(point(0.5, 1), poly2) == inside);

    // alternative version 2
    polygon poly3 = {
      point(2, 1),
      point(3, 1),
      point(4, 1),
      point(4, 0),
      point(0, 0),
      point(0, 2),
      point(1, 1),
      point(2, 1)
    };

    expect_true(point_in_polygon(point(-1, 1), poly3) == outside);
    expect_true(point_in_polygon(point(5, 1), poly3) == outside);
    expect_true(point_in_polygon(point(0.5, 1), poly3) == inside);

    // alternative version 3
    polygon poly4 = {
      point(4, 1),
      point(4, 0),
      point(0, 0),
      point(0, 2),
      point(1, 1),
      point(2, 1),
      point(3, 1),
      point(4, 1)
    };

    expect_true(point_in_polygon(point(-1, 1), poly4) == outside);
    expect_true(point_in_polygon(point(5, 1), poly4) == outside);
    expect_true(point_in_polygon(point(0.5, 1), poly4) == inside);
  }

  test_that("Zigzag 1") {
    polygon poly = {
      point(0, 2),
      point(1, 1),
      point(2, 1.5),
      point(3, 1),
      point(4, 1.5),
      point(5, 0),
      point(0, 0),
      point(0, 2)
    };

    expect_true(point_in_polygon(point(-1, 1), poly) == outside);
    expect_true(point_in_polygon(point(5, 1), poly) == outside);
    expect_true(point_in_polygon(point(0.5, 1), poly) == inside);
    expect_true(point_in_polygon(point(3, 1), poly) == undetermined);
  }

  test_that("Zigzag 2") {
    polygon poly = {
      point(0, 2),
      point(1, 1),
      point(2, 1.5),
      point(3, 1),
      point(4, 1.5),
      point(4, 3),
      point(0, 3),
      point(0, 2)
    };

    expect_true(point_in_polygon(point(-1, 1), poly) == outside);
    expect_true(point_in_polygon(point(5, 1), poly) == outside);
    expect_true(point_in_polygon(point(0.5, 1), poly) == outside);
    expect_true(point_in_polygon(point(1, 1.3), poly) == inside);
    expect_true(point_in_polygon(point(3, 1), poly) == undetermined);
  }

}

context("Polygon in polygon") {
  test_that("Basic relationships") {
    polygon p1 = {
      point(0, 0),
      point(0, 2),
      point(2, 2),
      point(2, 0),
      point(0, 0)
    };
    polygon p2 = {
      point(0.5, 0.5),
      point(0.5, 1.5),
      point(1.5, 1.5),
      point(1.5, 0.5),
      point(0.5, 0.5)
    };
    polygon p3 = {
      point(-1, -1),
      point(-1, 0),
      point(0, 0),
      point(0, -1),
      point(-1, -1)
    };
    polygon p4 = {
      point(-1, -1),
      point(-1, 1),
      point(1, 1),
      point(1, -1),
      point(-1, -1)
    };

    expect_true(polygon_in_polygon(p2, p1) == inside);
    expect_true(polygon_in_polygon(p1, p2) == outside);
    expect_true(polygon_in_polygon(p1, p3) == outside);
    expect_true(polygon_in_polygon(p3, p1) == outside);
    expect_true(polygon_in_polygon(p1, p4, false) == undetermined);
    expect_true(polygon_in_polygon(p4, p1, false) == undetermined);
  }

  test_that("Degenerate case") {
    polygon p1 = {
      point(0, 0),
      point(0, 2),
      point(2, 2),
      point(2, 0),
      point(0, 0)
    };
    expect_true(polygon_in_polygon(p1, p1) == undetermined);
  }
}


context("is_valid_ring()") {
  test_that("valid ring") {
    point p(0, 0);
    polygon poly;

    expect_false(is_valid_ring(poly));

    poly.push_back(p);
    expect_false(is_valid_ring(poly));

    poly.push_back(p);
    expect_false(is_valid_ring(poly));

    poly.push_back(p);
    expect_false(is_valid_ring(poly));

    poly.push_back(p);
    expect_false(is_valid_ring(poly));

    poly.push_back(point(1, 1));
    expect_true(is_valid_ring(poly));

    polygon poly2 = {
      point(0, 0),
      point(0, 2),
      point(2, 2),
      point(2, 0),
      point(0, 0)
    };
    expect_true(is_valid_ring(poly2));
  }
}
