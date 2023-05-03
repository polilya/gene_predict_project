/* A set of features based on DNA k-mer pairs */

#ifndef DNAKMERPAIR_H
#define DNAKMERPAIR_H

#include "SequenceFeatureSet.h"

class DNAKmerPairFS : public SequenceFeatureSet {
 public:
  string firstSeqName;
  string secondSeqName;
  int firstSeqID;
  int secondSeqID;
  int k;
  int* kmerPairIndices[3];

  DNAKmerPairFS(int k, int frame, string parameterization, string firstSeqName, string secondSeqName);
  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estseq);
  virtual void runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				vector< vector<weight_t> >& runningSumMinus);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores);
  int kmerToParameter(char* kmer, int k);
};

#endif
