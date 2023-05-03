#include "IntegerSequence.h"

IntegerSequence::IntegerSequence(int numIDs, int length) {
  
  numberOfSequences = 2 * numIDs;
  this->length = length;
  sequenceArray = new char*[numberOfSequences];
  for (int i=0; i<numberOfSequences; i++)
    sequenceArray[i] = new char[length];
}

IntegerSequence::~IntegerSequence() {
  for (int i=0; i<numberOfSequences; i++)
    delete[] sequenceArray[i];
  delete[] sequenceArray;
}
