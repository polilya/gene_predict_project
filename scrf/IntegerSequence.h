#ifndef INTEGERSEQUENCE_H
#define INTEGERSEQUENCE_H

#include "Sequence.h"

class IntegerSequence : public Sequence {
 public:
  char** sequenceArray;

  IntegerSequence(int numIDs, int length);
  ~IntegerSequence();
};

#endif
