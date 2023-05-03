#include "MaskingFeature.h"

MaskingFeature::MaskingFeature() {
  type = MASKING;
  frame = -1;

  valueToParam.resize(1);
  valueToWeight.resize(1);
  valueToGlobalParamID.resize(1);

  numParams = 1;
  valueToParam[0] = 0;
}

void MaskingFeature::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;
}

void MaskingFeature::scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) {
  //add our score if we see a lowercase (soft-masked) character
  for (pos_t j=0; j<alignment->length; j++) {
    switch (alignment->sequenceArray[0][j]) {
    case 'a':
    case 'c':
    case 'g':
    case 't':
      plusScores[j] += valueToWeight[0];
      minusScores[j] += valueToWeight[0];
    }
  }
}
