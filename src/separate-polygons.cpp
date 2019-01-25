// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <testthat.h>

#include <iostream>
#include <vector>
#include <set>

using namespace std;

#include "polygon.h"
#include "separate-polygons.h"

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
          wrap_around = true; // should never get here, due to choice of istart
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


in_polygon_type polygon_in_polygon(const polygon &query, const polygon &reference, bool fast) {
  int ins = 0, out = 0, undet = 0;

  for (int i = 0; i < query.size()-1; i++) {
    switch(point_in_polygon(query[i], reference)) {
    case inside:
      ins += 1;
      break;
    case outside:
      out += 1;
      break;
    default:
      undet += 1;
    }

    // shortcut for faster classification: if at least one
    // non-ambiguous point is found, we know whether we're inside
    // or outside
    if (fast && (ins > 0 || out > 0)) break;
  }

  if (ins > 0 && out == 0) {
    return inside;
  }

  if (out > 0 && ins == 0) {
    return outside;
  }

  return undetermined;
}


class polygon_hierarchy {
private:
  // for each polygon, contains a set of exterior polygons
  vector<set<int>> ext_polygons;
  vector<bool> active_polygons;

public:
  polygon_hierarchy(int n) {
    ext_polygons.resize(n);
    active_polygons.resize(n);

    // initially, all polygons are active
    for (auto it = active_polygons.begin(); it != active_polygons.end(); it++) {
      *it = true;
    }
  }

  void print() {
    for (int i = 0; i < ext_polygons.size(); i++) {
      cout << "polygon " << i << " (active = " << active_polygons[i] << ")" << endl;
      cout << "  enclosing: ";
      for (auto it = ext_polygons[i].begin(); it != ext_polygons[i].end(); it++) {
        cout << (*it) << " ";
      }
      cout << endl;
    }
  }

  void set_exterior(int poly, int exterior) {
    ext_polygons[poly].insert(exterior);
  }

  void remove(int poly) {
    for (auto it = ext_polygons.begin(); it != ext_polygons.end(); it++) {
      it->erase(poly);
    }
  }

  // returns the next top level polygon found
  int top_level_poly() {
    int i = 0;
    do {
      if (active_polygons[i] && ext_polygons[i].size() == 0) {
        active_polygons[i] = false;
        break;
      }
      i++;
    } while (i < ext_polygons.size());
    if (i == ext_polygons.size()) {
      // we have run out of top-level polygons, hence we're done
      i = -1;
    }

    return i;
  }

  // find all holes belonging to polygon, remove them and the parent
  // polygon from the hierarchy, and return
  set<int> collect_holes(int poly) {
    set<int> holes;

    int i = 0;
    do {
      if (active_polygons[i] &&
          ext_polygons[i].size() == 1 &&
          ext_polygons[i].count(poly) == 1) {
        holes.insert(i);
        active_polygons[i] = false;
      }
      i++;
    } while (i < ext_polygons.size());

    for (auto it = holes.begin(); it != holes.end(); it++) {
      remove(*it);
    }
    remove(poly);

    return holes;
  }
};

NumericMatrix polygon_as_matrix(polygon p, bool reverse = false) {
  int n = p.size();

  NumericMatrix m(n, 2);

  if (reverse) {
    for (int i = n; i > 0; i--) {
      m(n-i, 0) = p[i-1].x;
      m(n-i, 1) = p[i-1].y;
    }
  } else {
    for (int i = 0; i < n; i++) {
      m(i, 0) = p[i].x;
      m(i, 1) = p[i].y;
    }
  }

  return m;
}

// [[Rcpp::export]]
List separate_polygons(const NumericVector &x, const NumericVector &y, const IntegerVector &id) {
  List out; // final result

  int n = x.size();
  if (n == 0) {
    // set sf classes
    out.attr("class") = CharacterVector::create("XY", "MULTIPOLYGON", "sfg");
    return out;
  }
  if (y.size() != n || id.size() != n) {
    stop("Inputs x, y, and id must be of the same length.");
  }

  // create polygons from input data
  vector<polygon> polys;
  int cur_id = id[0];
  int cur_poly = 0;
  polys.push_back(polygon());
  for (int i = 0; i<n; i++) {
    if (id[i] != cur_id) {
      // complete current polygon and start new one
      polys.push_back(polygon());
      cur_id = id[i];
      cur_poly += 1;
    }
    polys[cur_poly].push_back(point(x[i], y[i]));
  }

  // close all polygons if necessary
  for (auto it = polys.begin(); it != polys.end(); it++) {
    if (!(it->front() == it->back())) {
      it->push_back(it->front());
    }
  }

  // set up polygon hierarchy
  polygon_hierarchy hi(polys.size());
  for (int i = 0; i < polys.size(); i++)
    for (int j = 0; j < polys.size(); j++ ) {
      if (i == j) continue;

      in_polygon_type result = polygon_in_polygon(polys[i], polys[j]);
      //cout << "polygon " << i << " is " << result << " of polygon " << j << endl;

      if (result == inside) {
        hi.set_exterior(i, j);
      }
      else if (result == undetermined){
        stop("Found polygons without undefined interior/exterior relationship.");
      }
    }

  int next_poly = hi.top_level_poly();
  while(next_poly >= 0) {
    List rings;
    rings.push_back(polygon_as_matrix(polys[next_poly]));
    // cout << "top-level polygon: " << next_poly << endl;

    set<int> holes = hi.collect_holes(next_poly);
    for (auto it = holes.begin(); it != holes.end(); it++) {
      // we reverse holes so they run in the same direction as outer polygons
      rings.push_back(polygon_as_matrix(polys[*it], true));
      //cout << "  hole: " << (*it) << endl;
    }
    out.push_back(rings);
    next_poly = hi.top_level_poly();
    //hi.print();
  }

  // set sf classes
  out.attr("class") = CharacterVector::create("XY", "MULTIPOLYGON", "sfg");
  return(out);
}

// testing code
/*** R
m <- matrix(c(0, 0, 0, 0, 0, 0,
              0, 1, 1, 1, 1, 0,
              0, 1, 2, 2, 1, 0,
              0, 1, 2, 0, 1, 0,
              0, 1, 1, 1, 1, 0,
              0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)

z <- isobands(1:6, 1:6, m, 0.5, 1.5)
mp1 <- separate_polygons(z[[1]]$x, z[[1]]$y, z[[1]]$id)

m <- matrix(c(0, 0, 0, 0, 0, 0,
              0, 2, 2, 2, 2, 0,
              0, 2, 0, 0, 2, 0,
              0, 2, 0, 0, 2, 0,
              0, 2, 2, 2, 2, 0,
              0, 0, 0, 0, 0, 0), 6, 6, byrow = TRUE)

z <- isobands(1:6, 1:6, m, 0.5, 1.5)
mp2 <- separate_polygons(z[[1]]$x, z[[1]]$y, z[[1]]$id)
*/
