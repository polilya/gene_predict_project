/* 
   A feature which gives a weight based on the previous 
   segment label and the current segment label 
*/

#ifndef TRANSITIONFEATURE_H
#define TRANSITIONFEATURE_H

#include "Globals.h"
#include "SequenceFeatureSet.h"
#include <vector>

class TransitionFeature {
 public:
  int typeID;
  stateid_t fromState;
  stateid_t toState;
  strand_t strand;
  weight_t weight;          // this weight is just for convenience...SemiCRF::transitionWeights holds the master weights
  int globalParamID;        // global paramter index of this transition's weight
  vector<char*> fromEnd;    /* sets of allowed characters at the end of the fromState
			       when making this transition */
  vector<char*> toStart;    /* sets of allowed characters at the beginning of the toState
			       when making this transition */
  vector<char*> disallowedFromEnd;
  vector<char*> disallowedToStart;

  bool allowed(AlignmentSequence* alignment, pos_t pos);
};

#endif
