// This file implements the 2D isoline and isoband algorithms described
// here: https://en.wikipedia.org/wiki/Marching_squares
// Line segments and polygons are not merged.
// Written by Claus O. Wilke

#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <vector>
#include <utility> // for pair
using namespace std;

// linear interpolation of boundary intersections
double interpolate(double x0, double x1, double z0, double z1, double value) {
  double d = (value - z0) / (z1 - z0);
  double x = x0 + d * (x1 - x0);
  return x;
}

// calculates the central value of a given cell
double central_value(int row, int col, const NumericMatrix &m) {
  return (m(row, col) + m(row, col + 1) + m(row + 1, col) + m(row + 1, col + 1))/4;
}

// functions that calculate the intersection points on the four corners of a cell
pair<double, double> intersect_top(int row, int col, double value, 
                                   const NumericVector &x, const NumericVector &y, const NumericMatrix &m) {
  return make_pair(
    interpolate(x[col], x[col+1], m(row, col), m(row, col+1), value),
    y[row]
  );
}

pair<double, double> intersect_bottom(int row, int col, double value, 
                                      const NumericVector &x, const NumericVector &y, const NumericMatrix &m) {
  return make_pair(
    interpolate(x[col], x[col+1], m(row+1, col), m(row+1, col+1), value),
    y[row+1]
  );
}

pair<double, double> intersect_left(int row, int col, double value, 
                                    const NumericVector &x, const NumericVector &y, const NumericMatrix &m) {
  return make_pair(
    x[col],
    interpolate(y[row], y[row+1], m(row, col), m(row+1, col), value)
  );
}

pair<double, double> intersect_right(int row, int col, double value, 
                                     const NumericVector &x, const NumericVector &y, const NumericMatrix &m) {
  return make_pair(
    x[col+1],
    interpolate(y[row], y[row+1], m(row, col+1), m(row+1, col+1), value)
  );
}

// functions that return one of the four corners of a cell
pair<double, double> top_left(int row, int col, const NumericVector &x, const NumericVector &y) {
  return make_pair(x[col], y[row]);
}

pair<double, double> top_right(int row, int col, const NumericVector &x, const NumericVector &y) {
  return make_pair(x[col + 1], y[row]);
}

pair<double, double> bottom_left(int row, int col, const NumericVector &x, const NumericVector &y) {
  return make_pair(x[col], y[row + 1]);
}

pair<double, double> bottom_right(int row, int col, const NumericVector &x, const NumericVector &y) {
  return make_pair(x[col + 1], y[row + 1]);
}


// [[Rcpp::export]]
DataFrame single_contour_lines(const NumericVector &x, const NumericVector &y, const NumericMatrix &m, double value) {
  // Basic marching squares implementation
  // Doesn't currently merge individual line segments, just writes out all of them individually
  // https://en.wikipedia.org/wiki/Marching_squares
  int nrow = m.nrow();
  int ncol = m.ncol();

  if (x.size() != ncol) {stop("Number of x coordinates must match number of columns in density matrix.");}
  if (y.size() != nrow) {stop("Number of y coordinates must match number of rows in density matrix.");}
    
  LogicalMatrix binarized(nrow, ncol, static_cast<LogicalVector>(m >= value).begin());
  IntegerMatrix cells(nrow - 1, ncol - 1);
  
  for (int r = 0; r < nrow-1; r++) {
    for (int c = 0; c < ncol-1; c++) {
      int index = 8*binarized(r, c) + 4*binarized(r, c + 1) + 2*binarized(r+1, c+1) + 1*binarized(r + 1, c);
      
      // two-segment saddles
      if (index == 5 && (central_value(r, c, m) >= value)) {
        index = 10;
      } else if (index == 10 && (central_value(r, c, m) >= value)) {
        index = 5;
      }
      
      cells(r, c) = index;
    }
  }
  
  // make line segments
  vector<double> x0, x1, y0, y1;  // vectors holding resulting line segments
  pair<double, double> p;    // point for temporary storage
  for (int r = 0; r < nrow-1; r++) {
    for (int c = 0; c < ncol-1; c++) {
      switch(cells(r, c)) {
      case 0: break;
      case 1:
        p = intersect_left(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_bottom(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 2:
        p = intersect_right(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_bottom(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 3:
        p = intersect_left(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_right(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 4:
        p = intersect_top(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_right(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 5:
        // case 2
        p = intersect_right(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_bottom(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        // case 7
        p = intersect_top(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_left(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 6:
        p = intersect_top(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_bottom(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 7:
        p = intersect_top(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_left(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 8:
        p = intersect_top(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_left(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 9:
        p = intersect_top(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_bottom(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 10:
        // case 1
        p = intersect_left(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_bottom(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        // case 4 
        p = intersect_top(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_right(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 11:
        p = intersect_top(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_right(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 12:
        p = intersect_left(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_right(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 13:
        p = intersect_right(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_bottom(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      case 14:
        p = intersect_left(r, c, value, x, y, m);
        x0.push_back(p.first); y0.push_back(p.second);
        p = intersect_bottom(r, c, value, x, y, m);
        x1.push_back(p.first); y1.push_back(p.second);
        break;
      default: break; // catch everything, just in case
      }
    }
  }
  
  return DataFrame::create(_["x0"] = x0, _["x1"] = x1, _["y0"] = y0, _["y1"] = y1, _["level"] = NumericVector(x0.size(), value));
}

// convenience function to append polygon points to output vectors
void push_point(pair<double, double> p, int cur_id, vector<double> &x_out, vector<double> &y_out, vector<int> &id) {
  x_out.push_back(p.first);
  y_out.push_back(p.second);
  id.push_back(cur_id);
}

// [[Rcpp::export]]
DataFrame single_contour_bands(const NumericVector &x, const NumericVector &y, const NumericMatrix &m, double vlo, double vhi) {
  int nrow = m.nrow();
  int ncol = m.ncol();
  
  if (x.size() != ncol) {stop("Number of x coordinates must match number of columns in density matrix.");}
  if (y.size() != nrow) {stop("Number of y coordinates must match number of rows in density matrix.");}
  
  IntegerVector v(nrow*ncol);
  IntegerVector::iterator iv = v.begin();
  for (NumericMatrix::const_iterator it = m.begin(); it != m.end(); it++) {
    *iv = (*it >= vlo && *it < vhi) + 2*(*it >= vhi);
    iv++;
  }
  
  IntegerMatrix ternarized(nrow, ncol, v.begin());
  IntegerMatrix cells(nrow - 1, ncol - 1);
  
  for (int r = 0; r < nrow-1; r++) {
    for (int c = 0; c < ncol-1; c++) {
      int index = 27*ternarized(r, c) + 9*ternarized(r, c + 1) + 3*ternarized(r + 1, c + 1) + ternarized(r + 1, c);
      
      cells(r, c) = index;
    }
  }
  
  
  // make polygons
  vector<double> x_out, y_out; vector<int> id;  // vectors holding resulting polygons
  pair<double, double> p;   // point for temporary storage
  int cur_id = 1;           // id counter for the polygons
  
  // all polygons are drawn clockwise for easy merging later on
  for (int r = 0; r < nrow-1; r++) {
    for (int c = 0; c < ncol-1; c++) {
      switch(cells(r, c)) {
      // doing cases out of order, sorted by type, is easier to keep track of

      // no contour
      case 0: break;
      case 80: break;

      // single triangle
      case 1: // 0001
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 3: // 0010
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 9: // 0100
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 27: // 1000
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 53: // 1222
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 71: // 2122
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 77: // 2212
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 79: // 2221
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
        
      // single trapezoid
      case 78: // 2220
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 74: // 2202
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 62: // 2022
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 26: // 0222
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 2: // 0002
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 6: // 0020
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 18: // 0200
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 54: // 2000
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
        
      // single rectangle
      case 4: // 0011
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 12: // 0110
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 36: // 1100
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 28: // 1001
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 76: // 2211
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 68: // 2112
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 44: // 1122
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 52: // 1221
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 72: // 2200
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 56: // 2002
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 8: // 0022
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 24: // 0220
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
        
      // single square
      case 40: // 1111
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
        
      // single pentagon
      case 49: // 1211
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 67: // 2111
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 41: // 1112
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 43: // 1121
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 31: // 1011
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 13: // 0111
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 39: // 1110
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 37: // 1101
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 45: // 1200
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 15: // 0120
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 5: // 0012
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 55: // 2001
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 35: // 1022
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 65: // 2102
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 75: // 2210
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 25: // 0221
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 29: // 1002
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 63: // 2100
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 21: // 0210
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 7: // 0021
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 51: // 1220
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 17: // 0122
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 59: // 2012
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 73: // 2201
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
   
      // single hexagon
      case 22: // 0211
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 66: // 2110
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 38: // 1102
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 34: // 1021
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 58: // 2011
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 14: // 0112
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 42: // 1120
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 46: // 1201
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 64: // 2101
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 16: // 0121
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 32: // 1012
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;
      case 48: // 1210
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
        push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
        push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
        push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
        cur_id++;
        break;

      // 8-sided saddle
      case 60: // 2020
        {
          double vc = central_value(r, c, m);
          if (vc < vlo) {
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else if (vc >= vhi) {
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 20: // 0202
        {
          double vc = central_value(r, c, m);
          if (vc < vlo) {
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else if (vc >= vhi) {
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;

      // 6-sided saddle
      case 10: // 0101
        {
          double vc = central_value(r, c, m);
          if (vc < vlo) {
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 30: // 1010
        {
          double vc = central_value(r, c, m);
          if (vc < vlo) {
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 70: // 2121
        {
          double vc = central_value(r, c, m);
          if (vc >= vhi) {
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 50: // 1212
        {
          double vc = central_value(r, c, m);
          if (vc >= vhi) {
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
        
      // 7-sided saddle
      case 69: // 2120
        {
          double vc = central_value(r, c, m);
          if (vc >= vhi) {
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 61: // 2021
        {
          double vc = central_value(r, c, m);
          if (vc >= vhi) {
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 47: // 1202
        {
          double vc = central_value(r, c, m);
          if (vc >= vhi) {
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 23: // 0212
        {
          double vc = central_value(r, c, m);
          if (vc >= vhi) {
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 11: // 0102
        {
          double vc = central_value(r, c, m);
          if (vc < vlo) {
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 19: // 0201
        {
          double vc = central_value(r, c, m);
          if (vc < vlo) {
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 33: // 1020
        {
          double vc = central_value(r, c, m);
          if (vc < vlo) {
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(top_left(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
      case 57: // 2010
        {
          double vc = central_value(r, c, m);
          if (vc < vlo) {
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            cur_id++;
          } else {
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            push_point(intersect_bottom(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_left(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vhi, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_top(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(intersect_right(r, c, vlo, x, y, m), cur_id, x_out, y_out, id);
            push_point(bottom_right(r, c, x, y), cur_id, x_out, y_out, id);
            cur_id++;
          }
        }
        break;
        
      default:
        // should never get here
        cout << "Unknown case: " << cells(r, c) << endl;
        break;
      }
    }
  }

  
  return DataFrame::create(_["x"] = x_out, _["y"] = y_out, _["id"] = id);
}

// The following R code will be automatically run after compilation.

/*** //R
library(ggplot2)

contour_lines <- function(m, level = NULL) {
  if (is.null(level)) {
    level <- scales::pretty_breaks(10)(m)
  }
  purrr::map_dfr(level, ~single_contour_lines(1:ncol(m), nrow(m):1, m, .x))
}

df <- contour_lines(volcano)
p <- ggplot(df, ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1, color = level)) + 
  geom_segment() +
  scale_color_viridis_c()

df1 <- single_contour_bands(1:ncol(volcano), nrow(volcano):1, volcano, 120, 140)
df2 <- single_contour_bands(1:ncol(volcano), nrow(volcano):1, volcano, 150, 152)
df3 <- contour_lines(volcano, level = c(120, 140, 150, 152))
ggplot(mapping = aes(x, y, group = id)) + 
  geom_polygon(data = df1, color = "lightblue", fill = "lightblue") +
  geom_polygon(data = df2, color = "tomato", fill = "tomato") +
  geom_segment(data = df3, aes(x = x0, y = y0, xend = x1, yend = y1), size = 0.2, inherit.aes = FALSE)
 
microbenchmark::microbenchmark(
  grDevices::contourLines(1:ncol(volcano), 1:nrow(volcano), volcano, levels = 120),
  single_contour_lines(1:ncol(volcano), 1:nrow(volcano), volcano, 120),
  single_contour_bands(1:ncol(volcano), 1:nrow(volcano), volcano, 120, 140)
)

*/
