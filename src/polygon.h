#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <unordered_map>

using namespace std;

// point in x-y space
struct point {
  double x, y; // x and y coordinates

  point(double x_in, double y_in) : x(x_in), y(y_in) {}
};

bool operator==(const point &p1, const point &p2);

typedef vector<point> polygon;

enum in_polygon_type {
  inside,       // point is inside a polygon
  outside,      // point is outside a polygon
  undetermined // point lies right on the boundary
};

#endif // POLYGON_H
