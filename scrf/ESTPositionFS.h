#ifndef ESTPOSITIONFS_H
#define ESTPOSITIONFS_H

#include "SequenceFeatureSet.h"

class ESTPositionFS : public SequenceFeatureSet {
 public:
  int k;

  ESTPositionFS();
  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estseq);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores);
};

#endif
