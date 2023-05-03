/* A feature set for features that depend only on the sequence and
   a given position, not any of the labels */

#ifndef SEQUENCEFEATURESET_H
#define SEQUENCEFEATURESET_H

#include "FeatureSet.h"
#include "AlignmentSequence.h"
#include "NumericSequence.h"
#include "ESTSequence.h"
#include <set>

class SequenceFeatureSet : public FeatureSet {
 public:
  int frame;         /* what frame can this SequenceFeatureSet score? -1 means all frames */
  vector< vector<weight_t> > runningSumPlus;
  vector< vector<weight_t> > runningSumMinus;
  sfsclass_t type;
  AlignmentSequence* alignment;
  NumericSequence* numericSeq;
  ESTSequence* estSeq;
  vector<string> stateTypeNames; /* names of state types associated with this SFS */
  
  virtual void setSequences(AlignmentSequence* alignment,
			    NumericSequence* numericSeq,
			    ESTSequence* estSeq) = 0;
  void positionExpectationHelper(pos_t pos, strand_t strand, weight_t prob, 
				 vector<weight_t>& counts);
  virtual void scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) = 0;
  void positionInitialGuessHelper(pos_t pos, strand_t strand, vector<weight_t>& denominators);
  void getDenominators(int index, set<int>& denominators);
};

#endif
