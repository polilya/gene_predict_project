#include <math.h>
#include "ESTPositionFS.h"

ESTPositionFS::ESTPositionFS() {
  this->type = ESTPOSITION;
  frame = -1;

  int numValues = EST_CHARS;
  valueToParam.resize(numValues);
  valueToWeight.resize(numValues);
  valueToGlobalParamID.resize(numValues);
  numParams = numValues;
  for (int i=0; i<numValues; i++)
    valueToParam[i] = i;
}

void ESTPositionFS::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;
}

void ESTPositionFS::scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) {
  for (pos_t j=0; j<alignment->length; j++) {
    plusScores[j] += valueToWeight[EST_TO_INDEX[estSeq->sequence[j]]];
    minusScores[j] += valueToWeight[EST_TO_INDEX[estSeq->sequence[j]]];
  }
}

