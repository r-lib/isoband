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
    double t = p1.x/(p1.x - p2.x);
    double yint = p1.y + t*(p2.y - p1.y);
    double xint = 0;

    if (yint < 0) { // actually need intersection with lower boundary
      t = p1.y/(p1.y - p2.y);
      xint = p1.x + t*(p2.x - p1.x);
      yint = 0;
    } else if (yint > 1) { // actually need intersection with upper boundary
      t = (1-p1.y)/(p2.y - p1.y);
      xint = p1.x + t*(p2.x - p1.x);
      yint = 1;
    }

    return point(xint, yint);
  }

  // p1 is to the right of box
  if (p1.x >= 1) {
    // intersection with right boundary
    double t = (1 - p1.x)/(p2.x - p1.x);
    double yint = p1.y + t*(p2.y - p1.y);
    double xint = 1;

    if (yint < 0) { // actually need intersection with lower boundary
      t = p1.y/(p1.y - p2.y);
      xint = p1.x + t*(p2.x - p1.x);
      yint = 0;
    } else if (yint > 1) { // actually need intersection with upper boundary
      t = (1-p1.y)/(p2.y - p1.y);
      xint = p1.x + t*(p2.x - p1.x);
      yint = 1;
    }

    return point(xint, yint);
  }

  // directly below
  if (p1.y <= 0) {
    // intersection with lower boundary
    double t = p1.y/(p1.y - p2.y);
    double xint = p1.x + t*(p2.x - p1.x);
    double yint = 0;
    return point(xint, yint);
  }

  // intersection with upper boundary
  double t = (1-p1.y)/(p2.y - p1.y);
  double xint = p1.x + t*(p2.x - p1.x);
  double yint = 1;
  return point(xint, yint);
}


// calculates the two intersection points ip1 and ip2 (if they exist) of a
// line segment from p1 to p2 and the unit box, assuming both p1 and p2
// are outside the box. if the assumption isn't true, results are not reliable.
// returns false if no intersection exists.
bool double_intersection(const point &p1, const point &p2, point &ip1, point &ip2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;

  if (dx == 0) {
    if (dy == 0) return false; // degenerate case

    // vertical line
    // trivial cases have been exlcuded by calling function, therefore this is easy
    ip1.x = p1.x;
    ip2.x = p1.x;
    if (p1.y >= 1) {
      ip1.y = 1;
      ip2.y = 0;
    } else {
      ip1.y = 0;
      ip2.y = 1;
    }
    return true;
  } else if (dy == 0) {
    // horizontal line
    // trivial cases have been exlcuded by calling function, therefore this is easy
    ip1.y = p1.y;
    ip2.y = p1.y;
    if (p1.x >= 1) {
      ip1.x = 1;
      ip2.x = 0;
    } else {
      ip1.x = 0;
      ip2.x = 1;
    }
    return true;
  } else {
    // in the general case, we need to calculate the intersection points with all four edges
    // we start at the top and go around clockwise, top, right, bottom, left
    double t1 = (1 - p1.y)/dy;
    double x1 = p1.x + t1*dx;
    double t2 = (1 - p1.x)/dx;
    double y2 = p1.y + t2*dy;
    double t3 = -1*p1.y/dy;
    double x3 = p1.x + t3*dx;
    double t4 = -1*p1.x/dx;
    double y4 = p1.y + t4*dy;

    int score = 8*(x1 > 0 && x1 <= 1) + 4*(y2 >= 0 && y2 < 1) + 2*(x3 >= 0 && x3 < 1) + (y4 > 0 && y4 <= 1);
    switch(score) {
    case 5: // 0101  line through left and right
      if (t4 < t2) {
        ip1 = point(0, y4);
        ip2 = point(1, y2);
      } else {
        ip1 = point(1, y2);
        ip2 = point(0, y4);
      }
      return true;
    case 10: // 1010  line through bottom and top
      if (t3 < t1) {
        ip1 = point(x3, 0);
        ip2 = point(x1, 1);
      } else {
        ip1 = point(x1, 1);
        ip2 = point(x3, 0);
      }
      return true;
    case 12: // 1100  line through top right
      if (t2 < t1) {
        ip1 = point(1, y2);
        ip2 = point(x1, 1);
      } else {
        ip1 = point(x1, 1);
        ip2 = point(1, y2);
      }
      return true;
    case 6: // 0110  line through bottom right
      if (t3 < t2) {
        ip1 = point(x3, 0);
        ip2 = point(1, y2);
      } else {
        ip1 = point(1, y2);
        ip2 = point(x3, 0);
      }
      return true;
    case 3: // 0011  line through bottom left
      if (t4 < t3) {
        ip1 = point(0, y4);
        ip2 = point(x3, 0);
      } else {
        ip1 = point(x3, 0);
        ip2 = point(0, y4);
      }
      return true;
    case 9: // 1001  line through top left
      if (t4 < t1) {
        ip1 = point(0, y4);
        ip2 = point(x1, 1);
      } else {
        ip1 = point(x1, 1);
        ip2 = point(0, y4);
      }
      return true;
    default:
      return false;
    }
  }
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

  // final case is double intersection in middle or no intersection at all
  bool crop = double_intersection(p1, p2, crop1, crop2);
  if (crop) return in_middle;
  return none;
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
  case in_middle:
    cout << "crop result: in middle " << crop1 << " " << crop2 << endl;
    break;
  default:
    cout << "crop result: not implemented" << endl;
  }
}


/*** R
test(c(-.3, .5), c(.5, 1.3))
*/
