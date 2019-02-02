#ifndef LABELS_H
#define LABELS_H

#include "polygon.h"

enum segment_crop_type {
  none,         // segment wasn't cropped
  complete,     // entire segment is gone
  at_beginning,    // beginning of segment is gone
  at_end,          // end of segment is gone
  in_middle        // middle of segment is gone
};

// crops the line segment running from p1 to p2 to a unit box
segment_crop_type crop_to_unit_box(const point &p1, const point &p2, point &crop1, point &crop2);

#endif // LABELS_H
