/* A sequence that the sCRF works with
   For convenience, the Sequence class actually consists
   of multiple sequences that can be of different types */

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "Globals.h"

class Sequence {
 public:
  pos_t length;
};

#endif
