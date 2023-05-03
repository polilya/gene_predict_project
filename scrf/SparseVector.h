#ifndef SPARSE_VECTOR
#define SPARSE_VECTOR

#include <iostream>
#include <vector>
#include "Globals.h"

#ifdef MULTI
#include <mpi.h>
#endif

class SparseVector {
 public:
  vector<int> indices;
  vector<weight_t> values;

  int size();
  void clear();
  void push_back(int index, weight_t value);
};

ostream& operator <<(ostream& os, const SparseVector& sv);

#endif
