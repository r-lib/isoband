// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <testthat.h>

#include <iostream>
#include <vector>
using namespace std;

#include "polygon.h"
#include "separate-polygons.h"

// eventually, move this to a polygon.cpp

ostream & operator<<(ostream &out, const point &p) {
  out << "(" << p.x << ", " << p.y << ")";
  return out;
}

ostream & operator<<(ostream &out, const in_polygon_type &t) {
  switch(t) {
  case inside:
    out << "inside";
    break;
  case outside:
    out << "outside";
    break;
  default:
    out << "undetermined";
  }
  return out;
}



/* Calculate the number of times a ray extending from point P to the right
 * intersects with the line segment defined by p0, p1. This number is
 * 0 or 1. However, -1 is returned if the point lies exactly on the segment,
 * so intersection in indetermined.
 */
int ray_intersections(point P, point p0, point p1) {
  checkUserInterrupt();
  //cout << P << " " << p0 << " " << p1 << endl;

  // simple cases
  if (p0.y < p1.y) {
    if ((P.y < p0.y) || (P.y > p1.y)) return 0;
  } else {
    if ((P.y > p0.y) || (P.y < p1.y)) return 0;
  }

  if ((P.x > p0.x) && (P.x > p1.x)) return 0;

  double dy = p1.y-p0.y;
  if (dy == 0) {
    if (P.y == p0.y) {
      // point is on the same y value, but does it lie inside the x interval?
      if ((P.x < p0.x) && (P.x < p1.x)) return 1;
      else return -1;
    }
    else return 0; // should never get here; handled by simple cases above
  }

  double t = (P.y - p0.y)/dy;
  double xint = p0.x + t*(p1.x - p0.x);
  //cout << "t = " << t << "; xint = " << xint << endl;
  if (xint < P.x) {
    return 0;
  }
  else if (xint == P.x) {
    return -1;
  }
  else return 1;
}

in_polygon_type point_in_polygon(const point &P, const polygon &poly) {
  int intersections = 0;
  int n = poly.size();
  int istart = 0;
  while (poly[istart].y == P.y) {
    // algorithm doesn't work if we start with a line segment that starts at P.y
    istart++;
    if (istart == n-1) {
      // degenerate polygon; one horizontal line
      // find min and max x and test if P.x is in
      // that interval or not
      double xmin = poly[0].x;
      double xmax = poly[0].x;
      for (int i = 1; i < n-1; i++) {
        if (poly[i].x < xmin) {
          xmin = poly[i].x;
        }
        if (poly[i].x > xmax) {
          xmax = poly[i].x;
        }
      }
      if (P.x >= xmin && P.x <= xmax) {
        return undetermined;
      } else {
        return outside;
      }
    }
  }

  int i = istart;
  do {
    int itr = ray_intersections(P, poly[i], poly[i+1]);
    //cout << i << " " << itr << endl;
    if (itr < 0) {
      // undetermined case, so we're done
      return undetermined;
    }

    if (itr > 0 && poly[i+1].y == P.y) {
      // special case, intersection with exact line endpoint
      bool from_above = poly[i].y > poly[i+1].y; // did we enter from above
      bool wrap_around = false;
      int j = i+1;
      do { // find next line segment where we move away from that point
        if (j == n-1) {
          j = 0;
        }
        if (j == istart) {
          wrap_around = true;
        }
        if (ray_intersections(P, poly[j], poly[j+1]) < 0) {
          // if the point lies exactly on any of these segments the case is undetermined
          return undetermined;
        }
        j++;
      } while (poly[j].y == poly[i+1].y);

      //cout << from_above << " " << i+1 << " " << j << " " << poly[i+1] << " " << poly[j] << endl;
      if ((!from_above && poly[j].y < poly[i+1].y) ||
          (from_above && poly[j].y > poly[i+1].y)) {
        // incorrect intersection
        //cout << "incorrect intersection" << endl;
        itr = 0;
      }
      i = j; // fast forward
      if (wrap_around || i == istart) {
        //cout << "have wrapped around during fast forward" << endl;
        //cout << "increment intersections (wa) at " << i << " " << itr << " " << intersections << endl;
        intersections += itr;
        break;
      }
      i--; // decrement by one because it'll be incremented again below
    }
    //cout << "increment intersections (el) at " << i << " " << itr << " " << intersections << endl;
    intersections += itr;
    i++;
    if (i == n-1) i = 0;
  } while(i != istart);

  if (intersections % 2 == 1) return inside;
  return outside;
}

// [[Rcpp::export]]
void separate_polygons() {
  polygon poly = {
    point(0, 0),
    point(0, 0)
  };
  in_polygon_type result;

  result = point_in_polygon(point(0, 1), poly);
  cout << "result: " << result << endl;
}

// testing code
/*** R
separate_polygons()
*/
