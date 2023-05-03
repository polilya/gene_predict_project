#include "Segment.h"

bool Segment::operator== (const Segment& s) {
  if (start != s.start || end != s.end || state != s.state)
    return false;
  return true;
}
