/* A type of state in the sCRF.  All states of this
   type share length distribution and the same FeatureSets. */

#ifndef STATETYPE_H
#define STATETYPE_H

#include "Globals.h"
#include "SequenceFeatureSet.h"
#include "LengthFS.h"
#include <string>
#include <vector>

class StateType {
 public:
  string name;
  LengthFS* lengthFeatures;
  vector<SequenceFeatureSet*> sequenceFeatureSets;
  vector<weight_t> plusScores;
  vector<weight_t> minusScores;

  StateType();
  ~StateType();
  void scorePositions();
  void freePositionScores();
  int findMatchingDNAKmerFS(int k, string parameterization, string seqName);
  int findMatchingDNAKmerPairFS(int k, string parameterization, 
				string firstSeqName, string secondSeqName);
};

#endif
