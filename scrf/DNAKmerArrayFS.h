/* A set of features based on DNA k-mers */

#ifndef DNAKMERARRAYFS_H
#define DNAKMERARRAYFS_H

#include "SequenceFeatureSet.h"

class DNAKmerArrayFS : public SequenceFeatureSet {
 public:
  string seqName;
  int seqID;
  int k;
  int length;
  int offset;
  int* kmerIndices[3];

  DNAKmerArrayFS(int k, int length, int offset, string parameterization, string seqName);
  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estseq);
  void runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
			vector< vector<weight_t> >& runningSumMinus);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores);
};

#endif
