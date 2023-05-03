/* A set of features based on numbers in an NumericSequence */

#ifndef NUMERICFS_H
#define NUMERICFS_H

#include "SequenceFeatureSet.h"

class NumericFS : public SequenceFeatureSet {
 public:
  int id;

  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq);
  virtual void runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				vector< vector<weight_t> >& runningSumMinus);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores);
};

#endif
