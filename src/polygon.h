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


typedef vector<point> polygon;

#endif // POLYGON_H
