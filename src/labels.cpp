#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
using namespace std;

#include "labels.h"

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
    return at_beginning;
  }

  if (p2_inside) {
    // simple case 2: crop at end
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

  cout << p1 << "; " << p2 << ": ";
  switch(result) {
  case none:
    cout << "crop result: none" << endl;
    break;
  case complete:
      cout << "crop result: complete" << endl;
    break;
  case at_beginning:
    cout << "crop result: at beginning" << endl;
    break;
  case at_end:
    cout << "crop result: at end" << endl;
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
test(c(-.2, -.3), c(.4, .5))

*/
