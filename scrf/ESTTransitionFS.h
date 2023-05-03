#ifndef ESTTRANSITIONFS_H
#define ESTTRANSITIONFS_H

#include "SequenceFeatureSet.h"

class ESTTransitionFS : public SequenceFeatureSet {
 public:
  int k;

  ESTTransitionFS();
  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estseq);
  virtual void runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				vector< vector<weight_t> >& runningSumMinus);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores);
};

#endif
