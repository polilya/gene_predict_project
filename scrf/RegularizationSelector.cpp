#include "RegularizationSelector.h"

RegularizationSelector::RegularizationSelector(Optimizer* optimizer, vector<weight_t>& initialGuessParams,
					       ofstream& trainingLog) {
  this->optimizer = optimizer;
  this->initialGuessParams = initialGuessParams;
  this->trainingLog = &trainingLog;

  minFunctionValue = INFINITY;
}

weight_t RegularizationSelector::ComputeFunction  (const vector<weight_t> &x) {
  optimizer->C = exp(x[0]);
  cerr << "Regularization selector trying C = " << optimizer->C
       << "(log C = " << log(optimizer->C) << ")" << endl;;
  (*trainingLog) << "Regularization selector trying C = " << optimizer->C 
		 << "(log C = " << log(optimizer->C) << ")" << endl;

  vector<weight_t> learningRate(initialGuessParams.size());
  for (int i=0; i<learningRate.size(); i++)
    learningRate[i] = max ((weight_t)1e-4, (weight_t)1e-3 * initialGuessParams[i]);
  
  vector<weight_t> params = initialGuessParams;
  optimizer->rProp(params, learningRate, 1.1, 0.5, true);

  weight_t functionValue = optimizer->minHoldoutObjective;
  if (functionValue < minFunctionValue) {
    minFunctionValue = functionValue;
    bestParams = optimizer->minHoldoutParams;
  }

  return functionValue;
}
