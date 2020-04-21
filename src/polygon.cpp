#include <iostream>
using namespace std;

#include "polygon.h"

ostream & operator<<(ostream &out, const point &p) {
  out << "(" << p.x << ", " << p.y << ")";
  return out;
}

bool operator==(const point &p1, const point &p2) {
  return (p1.x == p2.x) && (p1.y == p2.y);
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
