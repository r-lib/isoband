#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include <cstddef>  // for size_t
using namespace std;

#include "clip-lines.h"


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

// helper function for crop_lines(); checks whether a single point is inside the unit box
bool in_unit_box(const point &p) {
  if (p.x > 0 && p.x < 1 && p.y > 0 && p.y < 1) return true;
  return false;
}

// helper function for crop_lines()
void record_points(NumericVector &x_out, NumericVector &y_out, IntegerVector &id_out,
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


//' Clip lines to the outside of a box
//'
//' Clip lines to the outside of a box. The box is specified via midpoint, width,
//' height, and a rotation angle in radians. This is used to create space within
//' isolines for text labels or other annotations.
//'
//' @param x Numeric vector of x coordinates
//' @param y Numeric vector of y coordinates
//' @param id Integer vector of id numbers indicating which lines are connected
//' @param p_mid Numeric vector specifying x and y position of the box midpoint
//' @param width Box width
//' @param height Box height
//' @param theta Box angle, in radians
//' @export
// [[Rcpp::export]]
List clip_lines_impl(const NumericVector &x, const NumericVector &y, const IntegerVector &id,
                     const double p_mid_x, const double p_mid_y, const double width,
                     const double height, const double theta) {
  // output variables
  NumericVector x_out, y_out;
  IntegerVector id_out;

  // input checks
  if (x.size() != y.size()) {
    stop("Number of x and y coordinates must match.");
  }
  if (x.size() != id.size()) {
    stop("Number of x coordinates and id values must match.");
  }
  if (x.size() == 0) {
    // empty input, return empty output
    return List::create(_["x"] = x_out, _["y"] = y_out, _["id"] = id_out);
  }


  // set up transformation
  // lower left point of cropping rectangle
  point ll(p_mid_x - width*cos(theta)/2 + height*sin(theta)/2,
           p_mid_y - width*sin(theta)/2 - height*cos(theta)/2);
  // lower right point
  point lr(ll.x + width*cos(theta), ll.y + width*sin(theta));
  // upper right point
  point ul(ll.x - height*sin(theta), ll.y + height*cos(theta));

  //cout << "c(" << p_mid[0] << ", " << ll.x << ", " << lr.x << ", " << ul.x << ")/5, c(" <<
  //  p_mid[1] << ", " << ll.y << ", " << lr.y << ", " << ul.y << ")/5" << endl;

  unitbox_transformer t(ll, lr, ul);

  // crop
  int cur_id = id[0];
  int cur_id_out = 0; // first output id - 1
  point p1, p2, p1t, p2t;
  point crop1, crop2;
  p1 = point(x[0], y[0]);
  p1t = t.transform(p1);

  bool p1_recorded = in_unit_box(p1t); // record only if not in unit box, catches singlets
  bool p2_recorded = true; // when we first enter the loop, have only p1 unrecorded
  bool new_line_segment = true;

  size_t i = 1;
  while(i < x.size()) {
    if (cur_id != id[i]) {
      // id mismatch means we are starting a new line segment

      // first record any points that haven't been recorded yet. catches singlets
      record_points(x_out, y_out, id_out, p1, p2, cur_id_out,
                    p1_recorded, p2_recorded, new_line_segment);
      // now set up next line segment
      p1 = point(x[i], y[i]);
      p1t = t.transform(p1);
      cur_id = id[i];
      p1_recorded = in_unit_box(p1t); // record only if not in unit box, catches singlets
      new_line_segment = true;
      i++;
      continue;
    }
    p2 = point(x[i], y[i]);
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

  return List::create(_["x"] = x_out, _["y"] = y_out, _["id"] = id_out);
}


/*** R
# TODO:
# 1. regression tests:
#    - rotated clipping area for non-45 degree angle
# 2. clipping to multiple boxes at once

x <- c(0, 0, 1, 1, 0, 2, 3, 2.5, 2)
y <- c(0, 1, 1, 0, 0, 2, 2, 3, 2)
id <- c(1, 1, 1, 1, 1, 2, 2, 2, 2)
out <- clip_lines_impl(x, y, id, 1, 1, 1, 1, 0)
grid.newpage()
grid.polyline(x = out$x/5, y = out$y/5, id = out$id)
*/
