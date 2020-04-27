#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>
#include <testthat.h>

#include <iostream>
#include <vector>
#include <set>

using namespace std;

#include "polygon.h"
#include "separate-polygons.h"
#include "utils.h"

/* Calculate the number of times a ray extending from point P to the right
 * intersects with the line segment defined by p0, p1. This number is
 * 0 or 1. However, -1 is returned if the point lies exactly on the segment,
 * so intersection in indetermined.
 */
int ray_intersections(point P, point p0, point p1) {
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

  for (unsigned int i = 0; i < query.size()-1; i++) {
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
    for (unsigned int i = 0; i < ext_polygons.size(); i++) {
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
    unsigned int i = 0;
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

    unsigned int i = 0;
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


bool is_valid_ring(const polygon &poly) {
  if (poly.size() < 4) return false; // any polygon with fewer than four points is not a valid ring

  const point &p1 = poly.front();
  auto it = poly.begin();
  for (it++; it != poly.end(); it++) {
    if (!(p1 == *it)) {
      return true; // at least one point is different; we call it a valid ring
    }
  }

  return false; // degenerate polygon; we ignore it
}

SEXP polygon_as_matrix(polygon p, bool reverse = false) {
  int n = p.size();

  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, n, 2));
  double* m_p = REAL(m);

  if (reverse) {
    for (int i = n; i > 0; i--) {
      m_p[n-i] = p[i-1].x;
      m_p[n-i+n] = p[i-1].y;
    }
  } else {
    for (int i = 0; i < n; i++) {
      m_p[i] = p[i].x;
      m_p[i+n] = p[i].y;
    }
  }

  UNPROTECT(1);
  return m;
}

extern "C" SEXP separate_polygons(SEXP x, SEXP y, SEXP id) {

  BEGIN_CPP
  SEXP out; // final result
  SEXP cl = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(cl, 0, Rf_mkChar("XY"));
  SET_STRING_ELT(cl, 1, Rf_mkChar("MULTIPOLYGON"));
  SET_STRING_ELT(cl, 2, Rf_mkChar("sfg"));

  int n = Rf_length(x);
  if (n == 0) {
    // set sf classes
    out = PROTECT(Rf_allocVector(VECSXP, 0));
    Rf_classgets(out, cl);
    UNPROTECT(2);
    return out;
  }
  if (Rf_length(y) != n || Rf_length(id) != n) {
    Rf_error("Inputs x, y, and id must be of the same length.");
  }

  double* x_p = REAL(x);
  double* y_p = REAL(y);
  int* id_p = INTEGER(id);
  // create polygons from input data
  vector<polygon> polys;
  int cur_id = id_p[0];
  int cur_poly = 0;
  polys.push_back(polygon());
  for (int i = 0; i<n; i++) {
    if (id_p[i] != cur_id) {
      // complete current polygon and start new one
      polys.push_back(polygon());
      cur_id = id_p[i];
      cur_poly += 1;
    }
    polys[cur_poly].push_back(point(x_p[i], y_p[i]));
  }

  // close all polygons if necessary
  for (auto it = polys.begin(); it != polys.end(); it++) {
    if (!(it->front() == it->back())) {
      it->push_back(it->front());
    }
  }

  // set up polygon hierarchy
  polygon_hierarchy hi(polys.size());
  for (unsigned int i = 0; i < polys.size(); i++) {
    if (checkInterrupt())  {
      longjump_interrupt();
    }

    for (unsigned int j = 0; j < polys.size(); j++ ) {
      if (i == j) continue;

      in_polygon_type result = polygon_in_polygon(polys[i], polys[j]);
      //cout << "polygon " << i << " is " << result << " of polygon " << j << endl;

      if (result == inside) {
        hi.set_exterior(i, j);
      }
      else if (result == undetermined){
        Rf_error("Found polygons without undefined interior/exterior relationship.");
      }
    }
  }

  int next_poly = hi.top_level_poly();
  int i = 0;
  CollectorList all_rings;
  while(next_poly >= 0) {
    if (i % 1000 == 0 && checkInterrupt()) {
      longjump_interrupt();
    }
    i++;

    // for simplicity, we collect the rings even if the polygon
    // is not valid; we just keep track of this and ignore it at
    // the end; this reduces the risk of bugs
    bool valid_poly = is_valid_ring(polys[next_poly]);

    // collect the holes, if any
    set<int> holes = hi.collect_holes(next_poly);

    // record the polygon if valid
    if (valid_poly) {
      // collect all the rings belonging to this polygon
      SEXP rings = PROTECT(Rf_allocVector(VECSXP, holes.size() + 1));

      // collect the outer ring
      SET_VECTOR_ELT(rings, 0, polygon_as_matrix(polys[next_poly]));

      int k = 1;

      for (auto it = holes.begin(); it != holes.end(); it++) {
        if (is_valid_ring(polys[*it])) {
          // we reverse holes so they run in the same direction as outer polygons
          SET_VECTOR_ELT(rings, k, polygon_as_matrix(polys[*it], true));
          k++;
        }
      }

      // Shrink list to actual size because some holes may have been invalid
      rings = PROTECT(Rf_lengthgets(rings, k));
      all_rings.push_back(rings);
      UNPROTECT(2);
    }
    next_poly = hi.top_level_poly();
  }

  out = all_rings;
  Rf_classgets(out, cl);

  UNPROTECT(1);
  return(out);

  END_CPP
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
