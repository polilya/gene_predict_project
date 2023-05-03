/* A feature to score a masked position */

#ifndef MASKINGFEATURE_H
#define MASKINGFEATURE_H

#include "SequenceFeatureSet.h"

class MaskingFeature : public SequenceFeatureSet {
 public:
  MaskingFeature();
  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estseq);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores);
};

#endif
