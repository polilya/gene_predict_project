#ifndef NUMERICSEQUENCE_H
#define NUMERICSEQUENCE_H

#include "Sequence.h"
#include <vector>

class NumericSequence : public Sequence {
 public:
  float** sequenceArray;
  char** tags;
  int numberOfSequences;

  char** compressedSequenceArray;
  char** compressedTags;
  vector<unsigned long> compressedSequenceLengths;
  vector<unsigned long> compressedTagsLengths;

  void compressSequences();
  void uncompressSequences();
  void freeUncompressedSequences();
  void freeCompressedSequences();
  void allocateUncompressedSequences();


  NumericSequence(int numIDs, int length);
  ~NumericSequence();
};

#endif
