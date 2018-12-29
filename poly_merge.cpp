// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <vector>
#include <utility> // for pair
using namespace std;


enum point_type {
  grid,  // point on the original data grid
  hintersect_lo, // intersection with horizontal edge, low value
  hintersect_hi, // intersection with horizontal edge, high value
  vintersect_lo, // intersection with vertical edge, low value
  vintersect_hi  // intersection with vertical edge, high value
};

struct grid_point {
  int r, c; // row and column
  point_type type; // point type
  
  grid_point(double r_in = 0, double c_in = 0, point_type type_in = grid) : r(r_in), c(c_in), type(type_in) {}
  grid_point(const grid_point &p) : r(p.r), c(p.c), type(p.type) {}
};

struct point {
  double x, y; // x and y coordinates
  
  point(double x_in, double y_in) : x(x_in), y(y_in) {}
};

bool operator==(const grid_point &p1, const grid_point &p2) {
  return (p1.r == p2.r) && (p1.c == p2.c) && (p1.type == p2.type);
}

ostream & operator<<(ostream &out, const grid_point &p) {
  out << "(" << p.c << ", " << p.r << ", " << p.type << ") ";
  return out;
}

struct poly_connect {
  grid_point prev, next;
  bool empty; // has this connection been defined yet?
  bool collected; // has this connection been collected into a final polygon?
  
  poly_connect() : empty(true), collected(false) {};
};

ostream & operator<<(ostream &out, const poly_connect &pc) {
  if (pc.empty) {out << "empty";}
  else {out << "prev: " << pc.prev << "; next: " << pc.next << " ";}
  return out;
}

class poly_bands {
private:
  int nrow, ncol; // numbers of rows and columns
  const NumericVector &grid_x,&grid_y;
  const NumericMatrix &grid_z;
  double vlo, vhi; // low and high cutoff values
  grid_point tmp_poly[8]; // no elementary polygon has more than 8 vertices
  poly_connect tmp_poly_connect[8];
  int tmp_poly_size; // current number of elements in tmp_poly
  
  vector<vector<vector<poly_connect> > > polygons;
  
  // internal member functions
  void poly_start(int r, int c, point_type type) { // start a new elementary polygon
    tmp_poly[0].r = r;
    tmp_poly[0].c = c;
    tmp_poly[0].type = type;
    
    tmp_poly_size = 1;
  }
  
  void poly_add(int r, int c, point_type type) { // add point to elementary polygon
    tmp_poly[tmp_poly_size].r = r;
    tmp_poly[tmp_poly_size].c = c;
    tmp_poly[tmp_poly_size].type = type;
    
    tmp_poly_size++;
  }
  
  // one edge case I hadn't considered: two polygons that touch diagonally, so that
  // they overlap in exactly one point. this needs to be handled appropriately and cannot in the current framework.
  
  void poly_merge() { // merge current elementary polygon to prior polygons
    cout << "before merging:" << endl;
    
    // first, we figure out the right connections for current polygon
    for (int i = 0; i < tmp_poly_size; i++) {
      // create defined state in tmp_poly_connect[]
      // for each point, find previous and next point in polygon
      tmp_poly_connect[i].empty = false;
      tmp_poly_connect[i].collected = false;
      tmp_poly_connect[i].next = tmp_poly[(i+1<tmp_poly_size) ? i+1 : 0];
      tmp_poly_connect[i].prev = tmp_poly[(i-1>=0) ? i-1 : tmp_poly_size-1];
      
      cout << tmp_poly[i] << tmp_poly_connect[i] << endl;
      
      // now merge with existing polygons if needed
      const grid_point &p = tmp_poly[i];
      if (!polygons[p.r][p.c][p.type].empty) { // point has been used before, need to merge polygons
        if (tmp_poly_connect[i].next == polygons[p.r][p.c][p.type].prev) {
          if (tmp_poly_connect[i].prev == polygons[p.r][p.c][p.type].next) {
            // if both prev and next cancel, point can be deleted
            tmp_poly_connect[i].empty = true;
          } else {
            // otherwise, merge in "next" direction
            tmp_poly_connect[i].next = polygons[p.r][p.c][p.type].next;
            tmp_poly_connect[i].empty = false;
          }
        } else if (tmp_poly_connect[i].prev == polygons[p.r][p.c][p.type].next) {
          // opposite is true, merge in "prev" direction
          tmp_poly_connect[i].prev = polygons[p.r][p.c][p.type].prev;
          tmp_poly_connect[i].empty = false;
        } else {
          // should never get here
          cerr << "Something is wrong in polygon merging; vertices overlap but edges do not." << endl;
        }
      }
    }
    
    cout << "after merging:" << endl;
    
    // then we copy the connections into the polygon matrix
    for (int i = 0; i < tmp_poly_size; i++) {
      const grid_point &p = tmp_poly[i];
      polygons[p.r][p.c][p.type] = tmp_poly_connect[i];
      cout << p << tmp_poly_connect[i] << endl;
    }
    
    cout << "new grid:" << endl;
    print_polygons_state();
  }
  
  void print_polygons_state() {
    for (int r = 0; r < nrow; r++) {
      for (int c = 0; c < ncol; c++) {
        for (int type = 0; type < 5; type++) {
          if (!polygons[r][c][type].empty) {
            cout << c << " " << r << " " << type << ": " << polygons[r][c][type] << endl;
          }
        }
      }
    }
    cout << endl;
  }
  
  // linear interpolation of boundary intersections
  double interpolate(double x0, double x1, double z0, double z1, double value) {
    double d = (value - z0) / (z1 - z0);
    double x = x0 + d * (x1 - x0);
    return x;
  }
  
  // calculate output coordinates for a given grid point
  point calc_point_coords(const grid_point &p) {
    switch(p.type) {
    case grid:
      return point(grid_x[p.c], grid_y[p.r]);
    case hintersect_lo: // intersection with horizontal edge, low value
      return point(interpolate(grid_x[p.c], grid_x[p.c+1], grid_z(p.r, p.c), grid_z(p.r, p.c+1), vlo), grid_y[p.r]);
    case hintersect_hi: // intersection with horizontal edge, high value
      return point(interpolate(grid_x[p.c], grid_x[p.c+1], grid_z(p.r, p.c), grid_z(p.r, p.c+1), vhi), grid_y[p.r]);
    case vintersect_lo: // intersection with vertical edge, low value
      return point(grid_x[p.c], interpolate(grid_y[p.r], grid_y[p.r+1], grid_z(p.r, p.c), grid_z(p.r+1, p.c), vlo));
    case vintersect_hi: // intersection with vertical edge, high value
      return point(grid_x[p.c], interpolate(grid_y[p.r], grid_y[p.r+1], grid_z(p.r, p.c), grid_z(p.r+1, p.c), vhi));
    }
  }

public:
  poly_bands(const NumericVector &x, const NumericVector &y, const NumericMatrix &z, double value_low, double value_high) :
    grid_x(x), grid_y(y), grid_z(z), vlo(value_low), vhi(value_high)
  {
    nrow = grid_z.nrow();
    ncol = grid_z.ncol();
    
    if (grid_x.size() != ncol) {stop("Number of x coordinates must match number of columns in density matrix.");}
    if (grid_y.size() != nrow) {stop("Number of y coordinates must match number of rows in density matrix.");}
    
    // set up matrix of polygon points
    polygons.resize(nrow);
    for (auto it = polygons.begin(); it != polygons.end(); it++) {
      (*it).resize(ncol);
      for (auto jt = (*it).begin(); jt != (*it).end(); jt++) {
        (*jt).resize(5); // five different point types defined in point_type enum
      }
    }
  }
  
  void make_contour_bands() {
    // setup matrix of ternarized cell representations
    IntegerVector v(nrow*ncol);
    IntegerVector::iterator iv = v.begin();
    for (NumericMatrix::const_iterator it = grid_z.begin(); it != grid_z.end(); it++) {
      *iv = (*it >= vlo && *it < vhi) + 2*(*it >= vhi);
      iv++;
    }
    
    IntegerMatrix ternarized(nrow, ncol, v.begin());
    IntegerMatrix cells(nrow - 1, ncol - 1);
    
    for (int r = 0; r < nrow-1; r++) {
      for (int c = 0; c < ncol-1; c++) {
        int index = 27*ternarized(r, c) + 9*ternarized(r, c + 1) + 3*ternarized(r + 1, c + 1) + ternarized(r + 1, c);
        cells(r, c) = index;
        cout << index << " ";
      }
      cout << endl;
    }
    
    // all polygons must be drawn clockwise for proper merging
    for (int r = 0; r < nrow-1; r++) {
      for (int c = 0; c < ncol-1; c++) {
        switch(cells(r, c)) {
        // doing cases out of order, sorted by type, is easier to keep track of
        
        // no contour
        case 0: break;
        case 80: break;
          
          // single triangle
        case 1: // 0001
          poly_start(r, c, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 3: // 0010
          poly_start(r, c+1, vintersect_lo);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_lo);
          poly_merge();
          break;
        case 9: // 0100
          poly_start(r, c, hintersect_lo);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_lo);
          poly_merge();
          break;
        case 27: // 1000
          poly_start(r, c, vintersect_lo);
          poly_add(r, c, grid);
          poly_add(r, c, hintersect_lo);
          poly_merge();
          break;
          
          // single square
        case 40: // 1111
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        
        // single pentagon
        case 13: // 0111
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_merge();
          break;
        }
      }
    }
  }
  
  DataFrame collect_polygons() {
    // make polygons
    vector<double> x_out, y_out; vector<int> id;  // vectors holding resulting polygon paths
    int cur_id = 0;           // id counter for the polygon lines

    for (int r = 0; r < nrow; r++) {
      for (int c = 0; c < ncol; c++) {
        for (int type = 0; type < 5; type++) {
          if (!polygons[r][c][type].empty && !polygons[r][c][type].collected) { // skip any grid points that are empty or already collected
            // we have found a new polygon line; process it
            cur_id++;
            
            grid_point start = grid_point(r, c, static_cast<point_type>(type));
            grid_point cur = start;
            int i = 0;
            do {
              i++;
              point p = calc_point_coords(cur);
              x_out.push_back(p.x);
              y_out.push_back(p.y);
              id.push_back(cur_id);

              // record that we have processed this point              
              polygons[cur.r][cur.c][cur.type].collected = true;
              // advance to next point in polygon
              cur = polygons[cur.r][cur.c][cur.type].next;
            } while (!(cur == start) && i < 10000); // keep going until we reach the start point again
          }
        }
      }
    }
    return DataFrame::create(_["x"] = x_out, _["y"] = y_out, _["id"] = id);
  }
};


// [[Rcpp::export]]
DataFrame merged_contour_bands(const NumericVector &x, const NumericVector &y, const NumericMatrix &z, double value_low, double value_high) {
  poly_bands p(x, y, z, value_low, value_high);
  p.make_contour_bands();
  return p.collect_polygons();
}


/***R
m <- matrix(c(0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 0), 4, 4, byrow = TRUE)

m2 <- matrix(c(1, 1, 1, 1, 1, 1, 
              1, 1, 1, 1, 1, 1,
              1, 1, 0, 0, 1, 1,
              1, 1, 0, 0, 1, 1,
              1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1), 6, 6, byrow = TRUE)


df <- merged_contour_bands((1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1), m, .5, 1.5)
grid::grid.newpage()
grid::grid.path(df$x, df$y, df$id)
*/                         
                             