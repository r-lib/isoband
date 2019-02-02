#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
using namespace std;

#include "labels.h"


// calculates the intersection point of a line segment and
// the unit box, assuming p1 is outside and p2 is inside.
// if the assumption isn't true, results are not reliable.
point entry_intersection(const point &p1, const point &p2) {
  // p1 is to the left of box
  if (p1.x <= 0) {
    // intersection with left boundary
    double t = (p1.x - p2.x)/p1.x;
    double yint = p1.y + t*(p2.y - p1.y);
    double xint = 0;

    if (yint < 0) { // actually need intersection with lower boundary
      t = (p1.y - p2.y)/p1.y;
      xint = p1.x + t*(p2.x - p1.x);
      yint = 0;
    } else if (yint > 1) { // actually need intersection with upper boundary
      t = (p2.y - p1.y)/(1-p1.y);
      xint = p1.x + t*(p2.x - p1.x);
      yint = 1;
    }

    return point(xint, yint);
  }

  // p1 is to the right of box
  if (p1.x >= 1) {
    // intersection with right boundary
    double t = (p2.x - p1.x)/(1 - p1.x);
    double yint = p1.y + t*(p2.y - p1.y);
    double xint = 1;

    if (yint < 0) { // actually need intersection with lower boundary
      t = (p1.y - p2.y)/p1.y;
      xint = p1.x + t*(p2.x - p1.x);
      yint = 0;
    } else if (yint > 1) { // actually need intersection with upper boundary
      t = (p2.y - p1.y)/(1-p1.y);
      xint = p1.x + t*(p2.x - p1.x);
      yint = 1;
    }

    return point(xint, yint);
  }

  // directly below
  if (p1.y <= 0) {
    // intersection with lower boundary
    double t = (p1.y - p2.y)/p1.y;
    double xint = p1.x + t*(p2.x - p1.x);
    double yint = 0;
    return point(xint, yint);
  }

  // intersection with upper boundary
  double t = (p2.y - p1.y)/(1-p1.y);
  double xint = p1.x + t*(p2.x - p1.x);
  double yint = 1;
  return point(xint, yint);
}

segment_crop_type crop_to_unit_box(const point &p1, const point &p2, point &crop1, point &crop2) {
  // trivial case 1: line segment trivially outside box
  if ((p1.x <= 0 && p2.x <= 0) || (p1.x >= 1 && p2.x >= 1) ||
      (p1.y <= 0 && p2.y <= 0) || (p1.y >= 1 && p2.y >= 1)) {
    return none;
  }

  bool p1_inside = p1.x > 0 && p1.x < 1 && p1.y > 0 && p1.y < 1;
  bool p2_inside = p2.x > 0 && p2.x < 1 && p2.y > 0 && p2.y < 1;

  if (p1_inside) {
    // trivial case 2: line segment fully inside box
    if (p2_inside) {
      return complete;
    }

    // otherwise, simple case 1: crop at beginning
    crop1 = entry_intersection(p2, p1);
    return at_beginning;
  }

  if (p2_inside) {
    // simple case 2: crop at end
    crop1 = entry_intersection(p1, p2);
    return at_end;
  }


  return in_middle;
}

// [[Rcpp::export]]
void test(const NumericVector &p1, const NumericVector &p2) {
  point start(p1[0], p1[1]);
  point end(p2[0], p2[1]);
  point crop1, crop2;

  segment_crop_type result = crop_to_unit_box(start, end, crop1, crop2);

  cout << start << "; " << end << ": ";
  switch(result) {
  case none:
    cout << "crop result: none" << endl;
    break;
  case complete:
      cout << "crop result: complete" << endl;
    break;
  case at_beginning:
    cout << "crop result: at beginning " << crop1 << endl;
    break;
  case at_end:
    cout << "crop result: at end " << crop1 << endl;
    break;
  default:
    cout << "crop result: not implemented" << endl;
  }

}


/*** R
test(c(.2, .3), c(.4, .5))
test(c(-.2, .3), c(.4, -.5))
test(c(-.2, -.3), c(-.4, -.5))
test(c(.2, .3), c(1.5, .5))
test(c(1.5, .5), c(.5, .5))
test(c(.5, .5), c(1.5, 1.5))
*/
