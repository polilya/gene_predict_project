/* A set of features based on integers in an IntegerSequence */

#ifndef INTEGERFS_H
#define INTEGERFS_H

#include "SequenceFeatureSet.h"

class IntegerFS : public SequenceFeatureSet {
 public:
  int id;

  virtual void setSequences(AlignmentSequence* alignment, IntegerSequence* intSeq);
  virtual void runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				vector< vector<weight_t> >& runningSumMinus);
};

#endif
