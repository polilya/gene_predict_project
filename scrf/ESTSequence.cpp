#include <zlib.h>

#include "ESTSequence.h"

#define BUFFER_SIZE 1048576

ESTSequence::ESTSequence(string filename) {
  string estseq = "";
  char buffer[BUFFER_SIZE];

  gzFile estFile;
  estFile = gzopen (filename.c_str(), "rb");
  if (estFile == NULL) {
    cerr << "Could not open " << filename << endl;
    fatalError("Could not open estseq file");
  }
  while (! gzeof(estFile)) {
    int bytesRead = gzread (estFile, buffer, BUFFER_SIZE-1);
    buffer[bytesRead] = '\0';
    estseq += buffer;
  }
  gzclose (estFile);

  length = estseq.length();
  sequence = new char[length];
  memcpy (sequence, estseq.c_str(), length);
  reverseSequence = NULL;
  compressedSequence = NULL;

  compressSequence();
  freeUncompressedSequence();
}

ESTSequence::~ESTSequence() {
  freeUncompressedSequence();
  freeCompressedSequence();
}

void ESTSequence::freeUncompressedSequence() {
  if (sequence != NULL)
    delete[] sequence;
  sequence = NULL;
  if (reverseSequence != NULL)
    delete[] reverseSequence;
  reverseSequence = NULL;
}

void ESTSequence::freeCompressedSequence() {
  if (compressedSequence != NULL)
    delete[] compressedSequence;
  compressedSequence = NULL;
}

void ESTSequence::compressSequence() {
  if (compressedSequence != NULL)
    delete[] compressedSequence;

  unsigned long bufferSize = 2 * length;
  char* buffer = new char[bufferSize];

  compressedSequenceLength = bufferSize;
  compress2((Bytef*)buffer, &compressedSequenceLength, (Bytef*)sequence, length, 1);
  compressedSequence = new char[compressedSequenceLength];
  memcpy(compressedSequence, buffer, compressedSequenceLength);
  
  delete[] buffer;
}

void ESTSequence::uncompressSequence() {
  freeUncompressedSequence();

  sequence = new char[length];
  reverseSequence = new char[length];
  
  unsigned long uncompressedLength = length; 
  uncompress((Bytef*)sequence, &uncompressedLength, (Bytef*)compressedSequence, compressedSequenceLength);

  for (int i=0; i<length; i++)
    reverseSequence[i] = sequence[length - i - 1];
}
