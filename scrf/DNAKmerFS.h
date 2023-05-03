/* A set of features based on DNA k-mers */

#ifndef DNAKMERFS_H
#define DNAKMERFS_H

#include "SequenceFeatureSet.h"

class DNAKmerFS : public SequenceFeatureSet {
 public:
  string seqName;
  int seqID;
  int k;
  int* kmerIndices[3];

  DNAKmerFS(int k, int frame, string parameterization, string seqName);
  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estseq);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores);
};

#endif
