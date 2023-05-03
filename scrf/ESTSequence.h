#ifndef ESTSEQUENCE_H
#define ESTSEQUENCE_H

#include "Sequence.h"
#include <vector>

class ESTSequence : public Sequence {
 public:
  char* sequence;
  char* reverseSequence;
  char* compressedSequence;
  unsigned long length;
  unsigned long compressedSequenceLength;

  void compressSequence();
  void uncompressSequence();
  void freeCompressedSequence();
  void freeUncompressedSequence();
  ESTSequence(string filename);
  ~ESTSequence();
};

#endif
