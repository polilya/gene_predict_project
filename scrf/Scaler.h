#include "Optimizer.h"

class Scaler : public Minimizer {
 public:
  Optimizer* optimizer;
  vector<weight_t> initialGuessParams;

  Scaler(Optimizer* optimizer, vector<weight_t>& initialGuessParams);
  virtual weight_t ComputeFunction  (const vector<weight_t> &x);

};
