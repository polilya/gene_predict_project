/* A state in the sCRF */

#ifndef STATE_H
#define STATE_H

#include "Globals.h"
#include "FeatureSet.h"
#include "LengthFS.h"
#include <string>
#include <vector>

class TransitionFeature;  /* forward declaration */

class State {
 public:
  int typeID;
  string name;
  int phase;
  int initialFrame;               /* frame in which the state must begin.  -1 means no requirement. */
  strand_t strand;
  weight_t startWeight;           /* weight assigned to starting the parse in this state */
  int startParamID;               /* global parameter index of start weight */
  bool optimize;                  /* should we try to optimize performance for this state? */
  vector<int> transitionsTo;      /* list of all transition IDs to this state */
  vector<int> transitionsFrom;    /* list of all transition IDs from this state */
  vector<char*> disallowedKmers;  /* list of DNA K-mers that are not allowed to occur
				     in-frame in one of these states */
  ~State();
};

#endif
