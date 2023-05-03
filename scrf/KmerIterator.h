#ifndef KMERITERATOR_H
#define KMERITERATOR_H

#include "Globals.h"
#include <string>

class KmerIterator {
 public:
  KmerIterator(int k, string type);
  ~KmerIterator();
  void next();
  void init();

  int k;
  int* counter;
  char* kmer;
  int index;
  const char* alphabet;
  int alphabetSize;
  bool finished;
};

#endif
