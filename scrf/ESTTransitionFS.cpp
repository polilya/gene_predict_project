#include "ESTTransitionFS.h"
#include <math.h>

ESTTransitionFS::ESTTransitionFS() {
  this->type = ESTTRANSITION;
  frame = -1;

  int numValues = EST_CHARS * EST_CHARS;
  valueToParam.resize(numValues);
  valueToWeight.resize(numValues);
  valueToGlobalParamID.resize(numValues);
  numParams = numValues;
  for (int i=0; i<numValues; i++)
    valueToParam[i] = i;
}

void ESTTransitionFS::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;
}

void ESTTransitionFS::runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				 vector< vector<weight_t> >& runningSumMinus) {
}

void ESTTransitionFS::scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) {
}

