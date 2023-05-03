#ifndef DNAKMERPAIRARRAY_H
#define DNAKMERPAIRARRAY_H

#include "SequenceFeatureSet.h"

class DNAKmerPairArrayFS : public SequenceFeatureSet {
 public:
  string firstSeqName;
  string secondSeqName;
  int firstSeqID;
  int secondSeqID;
  int k;
  int length;
  int offset;
  int* kmerPairIndices[3];

  DNAKmerPairArrayFS(int k, int length, int offset, string parameterization, 
		     string firstSeqName, string secondSeqName);
  int kmerToParameter(char* kmer, int k);
  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estseq);
  void runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
			vector< vector<weight_t> >& runningSumMinus);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores);
};

#endif
