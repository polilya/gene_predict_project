#include "NumericFS.h"

void NumericFS::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;
}

void NumericFS::runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				 vector< vector<weight_t> >& runningSumMinus) {}

void NumericFS::scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) {}


