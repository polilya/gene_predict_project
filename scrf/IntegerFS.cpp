#include "IntegerFS.h"

void IntegerFS::setSequences(AlignmentSequence* alignment, IntegerSequence* intSeq) {
  this->alignment = alignment;
  this->intSeq = intSeq;
}

void IntegerFS::runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				 vector< vector<weight_t> >& runningSumMinus) {
  int plusRow = 2 * id;
  int minusRow = 2 * id + 1;

  vector<weight_t> ourRunningSumPlus(3, 0);
  vector<weight_t> ourRunningSumMinus(3, 0);

  for (pos_t j=0; j<intSeq->length; j++) {
    weight_t plusContribution = valueToWeight[intSeq->sequenceArray[plusRow][j]];
    weight_t minusContribution = valueToWeight[intSeq->sequenceArray[minusRow][alignment->length - j - 1]];

    if (frame == -1) {
      runningSumPlus[0][j] += plusContribution + ourRunningSumPlus[0];
      runningSumMinus[0][j] += minusContribution + ourRunningSumMinus[0];
      ourRunningSumPlus[0] += plusContribution;
      ourRunningSumMinus[0] += minusContribution;
    }
    else {
      fatalError("IntegerFS only supports frame of -1");
    }
  }
}
