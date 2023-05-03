#include "AlignmentBlock.h"
#include <cassert>
#include <zlib.h>

bool AlignmentBlock::operator<(const AlignmentBlock& other) const {
  return score < other.score;
} 

void AlignmentBlock::addSequence(int species, string& sequence) {
  this->species.push_back(species);
  assert(sequence.length() == length);
  sequences.push_back(new char[sequence.length()]);
  memcpy(sequences.back(), sequence.c_str(), length);
}

void AlignmentBlock::print(ostream& os) {
  os << "score: " << score << endl;
  os << "start: " << start << endl;
  os << species.size() << " species, " << sequences.size() << " sequences" << endl;
  for (int i=0; i<species.size(); i++) {
    os << species[i] << "\t" << sequences[i] << endl;
  }
}

void AlignmentBlock::compressBlock() {
  unsigned long bufferSize = 2 * length;
  char* buffer = new char[bufferSize];

  compressedSequenceLengths.resize(species.size());

  for (int i=0; i<species.size(); i++) {
    compressedSequenceLengths[i] = bufferSize;
    compress2((Bytef*)buffer, &(compressedSequenceLengths[i]), (Bytef*)sequences[i], length, 1);
    delete[] sequences[i];
    sequences[i] = new char[compressedSequenceLengths[i]];
    memcpy(sequences[i], buffer, compressedSequenceLengths[i]);
  }

  delete[] buffer;
}

void AlignmentBlock::uncompressBlock() {
  unsigned long uncompressedLength = length;
  unsigned long bufferSize = 2 * length;
  char* buffer = new char[bufferSize];

  for (int i = 0; i < species.size(); i++) {
    uncompress((Bytef*)buffer, &uncompressedLength, 
	       (Bytef*)sequences[i], compressedSequenceLengths[i]);
    delete[] sequences[i];
    sequences[i] = new char[length];
    memcpy(sequences[i], buffer, length);
  }

  delete[] buffer;
}
