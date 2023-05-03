#include "SparseVector.h"
#include <cassert>

int SparseVector::size() {
  return indices.size();
}

void SparseVector::clear() {
  vector<int> emptyIndices;
  vector<weight_t> emptyValues;

  indices.swap(emptyIndices);
  values.swap(emptyValues);
}

void SparseVector::push_back(int index, weight_t value) {
  assert(indices.size() == 0 || indices.back() < index);
  indices.push_back(index);
  values.push_back(value);
}

ostream& operator <<(ostream& os, const SparseVector& sv) {
  for (int i=0; i<sv.indices.size(); i++) {
    os << sv.indices[i] << ":" << sv.values[i];
    if (i < sv.indices.size() - 1)
      os << " ";
  }
  return os;
}
