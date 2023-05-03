#ifndef SEGMENT_H
#define SEGMENT_H

#include "Globals.h"

class Segment {
 public:
  stateid_t state;
  pos_t start;
  pos_t end;
  weight_t weight;

  bool operator == (const Segment& s);
};


#endif
