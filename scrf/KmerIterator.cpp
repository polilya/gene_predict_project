#include "KmerIterator.h"

KmerIterator::KmerIterator(int k, string type) {
  this->k = k;

  if (type == "DNA") {
    alphabet = INDEX_TO_DNA;
    alphabetSize = 7;
  }
  else
    fatalError("Unknown type for KmerIterator");

  counter = new int[k];
  kmer = new char[k+1];
  
  init();
}

void KmerIterator::init() {
  for (int i=0; i<k; i++) {
    counter[i] = 0;
    kmer[i] = alphabet[0];
  }
  kmer[k] = '\0';
  finished = false;
  index = 0;
}

void KmerIterator::next() {
  if (finished)
    return;

  /* increment counter, carry if needed */
  index++;
  counter[k-1]++;
  int i = k-1;
  while (i > 0 && counter[i] == alphabetSize) {
    counter[i] = 0;
    counter[i-1]++;
    i--;
  }
  if (counter[0] == alphabetSize) {
    finished = true;
    return;
  }

  /* set kmer string */
  for (i=0; i<k; i++)
    kmer[i] = alphabet[counter[i]];
  kmer[k] = '\0';
}

KmerIterator::~KmerIterator() {
  delete[] counter;
  delete[] kmer;
}
