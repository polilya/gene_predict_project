#include "FeatureSet.h"

void FeatureSet::setValueToWeight() {
  for (int i=0; i<valueToWeight.size(); i++) {
    valueToWeight[i] = paramToWeight[valueToParam[i]];
  }
}

void FeatureSet::setValueToGlobalParamID() {
  for (int i=0; i<valueToGlobalParamID.size(); i++) {
     valueToGlobalParamID[i] = globalParamID[valueToParam[i]];
  }
}
