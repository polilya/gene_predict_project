/* A type of transition in the sCRF.  All transitions of this type
   share the same SequenceFeatureSets that are used for scoring when
   the transition occurs. */

#ifndef TRANSITIONTYPE_H
#define TRANSITIONTYPE_H

#include <string>
#include <vector>
#include "Globals.h"
#include "SequenceFeatureSet.h"

class TransitionType {
 public:
  string name;
  vector<SequenceFeatureSet*> sequenceFeatureSets;

  ~TransitionType();
};



#endif
