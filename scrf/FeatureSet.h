/* A set of features used by the sCRF */

/* These are indicator features, that is, they take on
   value 1 if a particular part of the sequence matches a
   particular value, and 0 otherwise */
#ifndef FEATURESET_H
#define FEATURESET_H

#include "Globals.h"
#include "Sequence.h"
#include <string>
#include <vector>

class FeatureSet {
 public:
  vector<int> valueToParam;
  vector<weight_t> paramToWeight;
  vector<int> globalParamID;        /* keeps track of global index of each parameter */
  vector<weight_t> valueToWeight;   /* shortcut for value --> param --> weight */
  vector<int> valueToGlobalParamID; /* shortcut for value --> param --> global ID */
  string parameterization;          /* a name for how the parameter to weight mappings are set */
  int numParams;

  virtual void setValueToWeight();        /* fills in valueToWeight based on valueToParam and
					     paramToWeight */
  void setValueToGlobalParamID(); /* fills in valueToGlobalParamID based on valueToParam and
				     globalParamID */
};

#endif
