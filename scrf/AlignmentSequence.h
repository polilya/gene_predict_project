/* This class represents a multiple alignment of several species.
   For convenience, it also records the translation of the first sequence
   in each possible frame. */

#ifndef ALIGNMENTSEQUENCE_H
#define ALIGNMENTSEQUENCE_H

#include "Inline.h"
#include "Sequence.h"
#include <string>
#include <vector>
#include <set>

class AlignmentSequence : public Sequence {
 public:
  char** sequenceArray;
  int numberOfSequences;
  
  char** compressedSequenceArray;
  vector<unsigned long> compressedSequenceLengths;

  vector<string> species;
 
  int*** kmerIndices;
  int**** kmerPairIndices;
  int*** reverseKmerIndices;
  int**** reverseKmerPairIndices;
  int maxK;
  int maxSeq;
  int maxPairK;
  int maxPairSeq;

  AlignmentSequence(string alignFilename, set<string>& speciesToRead);
  ~AlignmentSequence();
  int getSpeciesID(string s);
  void reverseComplement(char* dna, char* rc);
  void write(string outFilename);
  void createReverseComplements();
  void freeReverseComplements();
  void maskLowercase();
  void compressSequences();
  void uncompressSequences();
  void freeUncompressedSequences();
  void freeCompressedSequences();
  void allocateUncompressedSequences();
  void freeKmerIndices();

  int getDNAKmerIndex(int seqNumber, pos_t pos, int k) {
    /*
    int index = 0;
    for (int i=0; i<k; i++)
      index += DNA_INDEX_COEFF[k - i - 1] * DNA_TO_INDEX[sequenceArray[seqNumber][pos + i]];
    return index;
    */

    return kmerIndices[k][seqNumber][pos];
  }

  int getDNAKmerPairIndex(int seqOneNumber, int seqTwoNumber, pos_t pos, int k) {
    /*
    int index = 0;
    for (int i=0; i<k; i++) {
      index += DNA_INDEX_COEFF[2*k - i - 1] * DNA_TO_INDEX[sequenceArray[seqOneNumber][pos + i]];
    }
    for (int i=0; i<k; i++) {
      index += DNA_INDEX_COEFF[k - i - 1] * DNA_TO_INDEX[sequenceArray[seqTwoNumber][pos + i]];
    }
    return index;    
    */

    return kmerPairIndices[k][seqOneNumber][seqTwoNumber][pos];
  }
};

#endif
