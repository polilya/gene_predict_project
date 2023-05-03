#include "Optimizer.h"

class RegularizationSelector : public Minimizer {
 public:
  Optimizer* optimizer;
  vector<weight_t> initialGuessParams;
  ofstream* trainingLog;

  vector<weight_t> bestParams;
  weight_t minFunctionValue;

  RegularizationSelector(Optimizer* optimizer, vector<weight_t>& initialGuessParams,
			 ofstream& trainingLog);
  virtual weight_t ComputeFunction  (const vector<weight_t> &x);
};
