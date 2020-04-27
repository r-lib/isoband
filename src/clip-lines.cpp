#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>

#include <iostream>
#include <cmath>
#include <cstddef>  // for size_t
using namespace std;

#include "clip-lines.h"
#include "utils.h"


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
    if (dy == 0) return false; // degenerate case, should never get here

    // vertical line
    // trivial cases have been excluded by calling function, therefore this is easy
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
    double t[4]; // linear parameters defining intersection points
    int b[4]; // integer to keep track of boundaries
    t[0] = (1 - p1.y)/dy; // top
    b[0] = 0;
    t[1] = (1 - p1.x)/dx; // right
    b[1] = 1;
    t[2] = -1*p1.y/dy; // bottom
    b[2] = 2;
    t[3] = -1*p1.x/dx; // left
    b[3] = 3;

    // now we need to sort the t values, we don't need
    // a complex algorithm here since it's so few cases
    for (int i = 1; i < 4; i++) { // find minimum
      if (t[0] > t[i]) {
        double temp = t[0];
        t[0] = t[i];
        t[i] = temp;

        int btemp = b[0];
        b[0] = b[i];
        b[i] = btemp;
      }
    }

    for (int i = 2; i < 4; i++) { // find next larger value
      if (t[1] > t[i]) {
        double temp = t[1];
        t[1] = t[i];
        t[i] = temp;

        int btemp = b[1];
        b[1] = b[i];
        b[i] = btemp;
      }
    }

    for (int i = 3; i < 4; i++) { // find next larger value
      if (t[2] > t[i]) {
        double temp = t[1];
        t[2] = t[i];
        t[i] = temp;

        int btemp = b[2];
        b[2] = b[i];
        b[i] = btemp;
      }
    }

    // t[1] and t[2] are the two inner-most intersections, which
    // define the two clipping points
    bool result = true;
    switch(b[1]) {
    case 0: // top
      ip1 = point(p1.x + t[1]*dx, 1);
      if (ip1.x < -1e-10 || ip1.x > 1+1e-10) result = false;
      break;
    case 1: // right
      ip1 = point(1, p1.y + t[1]*dy);
      if (ip1.y < -1e-10 || ip1.y > 1+1e-10) result = false;
      break;
    case 2: // bottom
      ip1 = point(p1.x + t[1]*dx, 0);
      if (ip1.x < -1e-10 || ip1.x > 1+1e-10) result = false;
      break;
    case 3: // left
      ip1 = point(0, p1.y + t[1]*dy);
      if (ip1.y < -1e-10 || ip1.y > 1+1e-10) result = false;
      break;
    default: // should never go here
      result = false;
    }

    switch(b[2]) {
    case 0: // top
      ip2 = point(p1.x + t[2]*dx, 1);
      if (ip2.x < -1e-10 || ip2.x > 1+1e-10) result = false;
      break;
    case 1: // right
      ip2 = point(1, p1.y + t[2]*dy);
      if (ip2.y < -1e-10 || ip2.y > 1+1e-10) result = false;
      break;
    case 2: // bottom
      ip2 = point(p1.x + t[2]*dx, 0);
      if (ip2.x < -1e-10 || ip2.x > 1+1e-10) result = false;
      break;
    case 3: // left
      ip2 = point(0, p1.y + t[2]*dy);
      if (ip2.y < -1e-10 || ip2.y > 1+1e-10) result = false;
      break;
    default: // should never go here
      result = false;
    }
    return result;
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

// helper function for crop_lines(); checks whether a single point is inside the unit box
bool in_unit_box(const point &p) {
  if (p.x > 0 && p.x < 1 && p.y > 0 && p.y < 1) return true;
  return false;
}

// helper function for crop_lines()
void record_points(vector<double> &x_out, vector<double> &y_out, vector<int> &id_out,
                   const point &p1, const point &p2, int &cur_id_out,
                   bool &p1_recorded, bool &p2_recorded, bool &new_line_segment) {
  if (new_line_segment) {
    // start a new line segment, but defer if nothing to record
    if (!p1_recorded || !p2_recorded) {
      cur_id_out++;
      new_line_segment = false;
    }
  }

  if (!p1_recorded) {
    x_out.push_back(p1.x);
    y_out.push_back(p1.y);
    id_out.push_back(cur_id_out);
    p1_recorded = true;
  }

  if (!p2_recorded) {
    x_out.push_back(p2.x);
    y_out.push_back(p2.y);
    id_out.push_back(cur_id_out);
    p2_recorded = true;
  }
}


// Clip lines to the outside of a box
//
// Clip lines to the outside of a box. The box is specified via midpoint, width,
// height, and a rotation angle in radians. This is used to create space within
// isolines for text labels or other annotations.
//
// @param x Numeric vector of x coordinates
// @param y Numeric vector of y coordinates
// @param id Integer vector of id numbers indicating which lines are connected
// @param p_mid_x,p_mid_y Numeric values specifying the x and y position of the box midpoint
// @param width Box width
// @param height Box height
// @param theta Box angle, in radians
// @param asp Aspect ratio (width/height) of the target canvas. This is used to convert widths
//  to heights and vice versa for rotated boxes
// @export
extern "C" SEXP clip_lines_impl(SEXP x, SEXP y, SEXP id, SEXP _p_mid_x, SEXP _p_mid_y,
                     SEXP _width, SEXP _height, SEXP _theta, SEXP _asp) {

  BEGIN_CPP
  // input
  int n = Rf_length(x);
  double* x_p = REAL(x);
  double* y_p = REAL(y);
  int* id_p = INTEGER(id);
  double p_mid_x = REAL(_p_mid_x)[0];
  double p_mid_y = REAL(_p_mid_y)[0];
  double width = REAL(_width)[0];
  double height = REAL(_height)[0];
  double theta = REAL(_theta)[0];
  double asp = REAL(_asp)[0];

  // output variable
  SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("x"));
  SET_STRING_ELT(names, 1, Rf_mkChar("y"));
  SET_STRING_ELT(names, 2, Rf_mkChar("id"));
  Rf_setAttrib(res, Rf_install("names"), names);

  vector<double> x_out, y_out;
  vector<int> id_out;

  // input checks
  if (n != Rf_length(y)) {
    Rf_error("Number of x and y coordinates must match.");
  }
  if (n != Rf_length(id)) {
    Rf_error("Number of x coordinates and id values must match.");
  }
  if (n == 0) {
    // empty input, return empty output
    SET_VECTOR_ELT(res, 0, Rf_allocVector(REALSXP, 0));
    SET_VECTOR_ELT(res, 1, Rf_allocVector(REALSXP, 0));
    SET_VECTOR_ELT(res, 2, Rf_allocVector(INTSXP, 0));
    UNPROTECT(2);
    return res;
  }

  // set up transformation
  // lower left point of cropping rectangle
  point ll(p_mid_x - width*cos(theta)/2 + (height/asp)*sin(theta)/2,
           p_mid_y - asp*width*sin(theta)/2 - height*cos(theta)/2);
  // lower right point
  point lr(ll.x + width*cos(theta), ll.y + asp*width*sin(theta));
  // upper left point
  point ul(ll.x - (height/asp)*sin(theta), ll.y + height*cos(theta));

  unitbox_transformer t(ll, lr, ul);

  // crop
  int cur_id = id_p[0];
  int cur_id_out = 0; // first output id - 1
  point p1, p2, p1t, p2t;
  point crop1, crop2;
  p1 = point(x_p[0], y_p[0]);
  p1t = t.transform(p1);

  bool p1_recorded = in_unit_box(p1t); // record only if not in unit box, catches singlets
  bool p2_recorded = true; // when we first enter the loop, have only p1 unrecorded
  bool new_line_segment = true;

  int i = 1;
  while(i < n) {
    if (cur_id != id_p[i]) {
      // id mismatch means we are starting a new line segment

      // first record any points that haven't been recorded yet. catches singlets
      record_points(x_out, y_out, id_out, p1, p2, cur_id_out,
                    p1_recorded, p2_recorded, new_line_segment);
      // now set up next line segment
      p1 = point(x_p[i], y_p[i]);
      p1t = t.transform(p1);
      cur_id = id_p[i];
      p1_recorded = in_unit_box(p1t); // record only if not in unit box, catches singlets
      new_line_segment = true;
      i++;
      continue;
    }
    p2 = point(x_p[i], y_p[i]);
    p2t = t.transform(p2);
    p2_recorded = false;
    segment_crop_type result = crop_to_unit_box(p1t, p2t, crop1, crop2);
    switch(result) {
    case complete:
      // skip recording for this line segment
      p1_recorded = true;
      p2_recorded = true;
      // start new line segment with next point
      new_line_segment = true;
      break;
    case at_beginning:
      p1t = crop1;
      p1 = t.inv_transform(p1t);
      p1_recorded = false;
      new_line_segment = true;
      break;
    case at_end:
      p2_recorded = false;
      record_points(x_out, y_out, id_out, p1, t.inv_transform(crop1), cur_id_out,
                    p1_recorded, p2_recorded, new_line_segment);
      new_line_segment = true;
      break;
    case in_middle:
      p2_recorded = false;
      record_points(x_out, y_out, id_out, p1, t.inv_transform(crop1), cur_id_out,
                    p1_recorded, p2_recorded, new_line_segment);
      p1t = crop2;
      p1 = t.inv_transform(p1t);
      p1_recorded = false;
      p2_recorded = false;
      new_line_segment = true;
      break;
    default:   // nothing to be done, record and move on
      break;
    }

    record_points(x_out, y_out, id_out, p1, p2, cur_id_out,
                  p1_recorded, p2_recorded, new_line_segment);
    p1 = p2;
    p1t = p2t;
    i++;
  }
  // record any remaining points; catches singlets
  record_points(x_out, y_out, id_out, p1, p2, cur_id_out,
                p1_recorded, p2_recorded, new_line_segment);

  int final_size = x_out.size();
  SEXP x_final = SET_VECTOR_ELT(res, 0, Rf_allocVector(REALSXP, final_size));
  double* x_final_p = REAL(x_final);
  SEXP y_final = SET_VECTOR_ELT(res, 1, Rf_allocVector(REALSXP, final_size));
  double* y_final_p = REAL(y_final);
  SEXP id_final = SET_VECTOR_ELT(res, 2, Rf_allocVector(INTSXP, final_size));
  int* id_final_p = INTEGER(id_final);

  for (int i = 0; i < final_size; ++i) {
    x_final_p[i] = x_out[i];
    y_final_p[i] = y_out[i];
    id_final_p[i] = id_out[i];
  }

  UNPROTECT(2);
  return res;

  END_CPP
}


/*** R
x <- c(0, 0, 1, 1, 0, 2, 3, 2.5, 2)
y <- c(0, 1, 1, 0, 0, 2, 2, 3, 2)
id <- c(1, 1, 1, 1, 1, 2, 2, 2, 2)
out <- clip_lines_impl(x, y, id, 1.5, 1.5, 2.5, 1, pi/4)
grid.newpage()
grid.polyline(x = out$x/5, y = out$y/5, id = out$id)
*/
