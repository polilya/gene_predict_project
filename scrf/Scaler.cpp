#include "Scaler.h"

Scaler::Scaler(Optimizer* optimizer, vector<weight_t>& initialGuessParams) {
  this->optimizer = optimizer;
  this->initialGuessParams = initialGuessParams;
}

weight_t Scaler::ComputeFunction  (const vector<weight_t> &x) {
  vector<weight_t> scaledParams = initialGuessParams;

  for (int i=0; i<optimizer->scrf->transitionTypes.size(); i++) {
    for (int k=0; k<optimizer->scrf->transitionTypes[i].sequenceFeatureSets.size(); k++) {
      SequenceFeatureSet& sfs = *(optimizer->scrf->transitionTypes[i].sequenceFeatureSets[k]);
      if (sfs.type == SVM) {
	for (int j=0; j<sfs.paramToWeight.size(); j++) {
	  scaledParams[sfs.globalParamID[j]] *= x[0];
	}
	if (optimizer->scrf->transitionTypes[i].name == "DonorSpliceGC") {
	  //take into account the fact that roughly 1% of donor sites are GC
	  for (int j=0; j<sfs.paramToWeight.size(); j++) {
	    scaledParams[sfs.globalParamID[j]] += log(0.01);
	  }
	}
      }
    }
  }

  return optimizer->ComputeFunction(scaledParams);
}
