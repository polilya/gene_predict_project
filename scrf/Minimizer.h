#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <vector>
#include <math.h>
#include <iostream>
#include <cassert>
#include "Globals.h"

using namespace std;

weight_t sign(weight_t x);
weight_t DotProduct (const vector<weight_t> &x, const vector<weight_t> &y);
weight_t Norm (const vector<weight_t> &x);
vector<weight_t>& Scale(vector<weight_t>& x, const vector<weight_t> &scale);
vector<weight_t>& Unscale(vector<weight_t>& x, const vector<weight_t> &scale);
vector<weight_t> operator+(const vector<weight_t> &x, const vector<weight_t> &y);
vector<weight_t> operator-(const vector<weight_t> &x, const vector<weight_t> &y);
vector<weight_t> &operator+= (vector<weight_t> &x, weight_t c);
vector<weight_t> operator*(const vector<weight_t> &x, weight_t c);

class Minimizer {
 public:
  virtual bool Report (const vector<weight_t> &theta, int iteration, weight_t objective, weight_t step_length);
  virtual void ComputeGradient (vector<weight_t> &g, const vector<weight_t> &x);
  virtual weight_t ComputeFunction  (const vector<weight_t> &x) = 0;
  virtual weight_t ComputeFunctionAndGradient (vector<weight_t> &g, const vector<weight_t> &x);
  Minimizer();
  virtual ~Minimizer(){}

  void rProp(vector<weight_t>& x, 
	     vector<weight_t>& learningRate,
	     weight_t etaPlus = 1.2, 
	     weight_t etaMinus = 0.5,
	     bool surge = false);

  weight_t GoldenSection(vector<weight_t>& x, vector<weight_t>& d, 
			 weight_t left, weight_t mid, weight_t right);
};

struct RpropState {
  vector<weight_t> x;
  weight_t functionValue;
  vector<weight_t> gradient;
  vector<weight_t> learningRates;
  vector<weight_t> metaLearningRates;
  vector<bool> holds;
};

#endif
