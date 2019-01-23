#ifndef SEPARATE_POLYGONS_H
#define SEPARATE_POLYGONS_H

#include "polygon.h"


/* Calculate the number of times a ray extending from point P to the right
 * intersects with the line segment defined by p0, p1. This number is
 * 0 or 1. However, -1 is returned if the point lies exactly on the segment,
 * so intersection in indetermined.
 */
int ray_intersections(point P, point p0, point p1);

/* Test whether a point lies inside a polygon or not. Can return one of
 * three values, inside, outside, or undetermined.
 */
in_polygon_type point_in_polygon(const point &P, const polygon &poly);

#endif
