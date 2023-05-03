#include <zlib.h>

#include "NumericSequence.h"

NumericSequence::NumericSequence(int numIDs, int length) {
  
  numberOfSequences = 2 * numIDs;
  this->length = length;

  compressedSequenceArray = NULL;
  compressedTags = NULL;
  compressedSequenceLengths.resize(numberOfSequences);
  compressedTagsLengths.resize(numberOfSequences);

  allocateUncompressedSequences();
}

NumericSequence::~NumericSequence() {
  freeUncompressedSequences();
  freeCompressedSequences();
}

void NumericSequence::compressSequences() {
  freeCompressedSequences();

  unsigned long bufferSize = 2 * sizeof(float) * length;
  char* buffer = new char[bufferSize];

  compressedSequenceArray = new char*[numberOfSequences];
  for (int i=0; i<numberOfSequences; i++) {
    compressedSequenceLengths[i] = bufferSize;
    compress2((Bytef*)buffer, &(compressedSequenceLengths[i]), (Bytef*)sequenceArray[i], sizeof(float)*length, 1);
    compressedSequenceArray[i] = new char[compressedSequenceLengths[i]];
    memcpy(compressedSequenceArray[i], buffer, compressedSequenceLengths[i]);
  }

  compressedTags = new char*[numberOfSequences];
  for (int i=0; i<numberOfSequences; i++) {
    compressedTagsLengths[i] = bufferSize;
    compress2((Bytef*)buffer, &(compressedTagsLengths[i]), (Bytef*)tags[i], length, 1);
    compressedTags[i] = new char[compressedTagsLengths[i]];
    memcpy(compressedTags[i], buffer, compressedTagsLengths[i]);
  }

  delete[] buffer;
}

void NumericSequence::uncompressSequences() {
  freeUncompressedSequences();
  allocateUncompressedSequences();
  
  unsigned long uncompressedLength = sizeof(float) * length; 

  for (int i=0; i<numberOfSequences; i++)
    uncompress((Bytef*)sequenceArray[i], &uncompressedLength, (Bytef*)compressedSequenceArray[i], compressedSequenceLengths[i]);

  uncompressedLength = length;

  for (int i=0; i<numberOfSequences; i++)
    uncompress((Bytef*)tags[i], &uncompressedLength, (Bytef*)compressedTags[i], compressedTagsLengths[i]);
}

void NumericSequence::freeUncompressedSequences() {
  if (sequenceArray != NULL) {
    for (int i=0; i<numberOfSequences; i++)
      delete[] sequenceArray[i];
    delete[] sequenceArray;
    sequenceArray = NULL;
  }

  if (tags != NULL) {
    for (int i=0; i<numberOfSequences; i++)
      delete[] tags[i];
    delete[] tags;
    tags = NULL;
  }
}

void NumericSequence::freeCompressedSequences() {
  if (compressedSequenceArray != NULL) {
    for (int i=0; i<numberOfSequences; i++)
      delete[] compressedSequenceArray[i];
    delete[] compressedSequenceArray;
    compressedSequenceArray = NULL;
  }

  if (compressedTags != NULL) {
    for (int i=0; i<numberOfSequences; i++)
      delete[] compressedTags[i];
    delete[] compressedTags;
  }
}

void NumericSequence::allocateUncompressedSequences() {
  sequenceArray = new float*[numberOfSequences];
  for (int i=0; i<numberOfSequences; i++)
    sequenceArray[i] = new float[length];

  tags = new char*[numberOfSequences];
  for (int i=0; i<numberOfSequences; i++)
    tags[i] = new char[length];
}
