#include "Minimizer.h"
#include <iomanip>
#include <cassert>

Minimizer::Minimizer() {
}

bool Minimizer::Report (const vector<weight_t> &theta, int iteration, weight_t objective, weight_t step_length) {
  cerr << "Report() called but was not overriden by a derived class";
  exit(0); 
}

void Minimizer::ComputeGradient (vector<weight_t> &g, const vector<weight_t> &x) {
  cerr << "ComputeGradient() called but was not overriden by a derived class";
  exit(0);
}

weight_t Minimizer::ComputeFunctionAndGradient (vector<weight_t> &g, const vector<weight_t> &x) {
  ComputeGradient(g, x);
  return ComputeFunction(x);
}

weight_t Minimizer::GoldenSection(vector<weight_t>& x, vector<weight_t>& d, 
				  weight_t left, weight_t mid, weight_t right) {
  const weight_t FINAL_BRACKET_SIZE = 0.1;
  const int MAX_ITERATIONS = 30;  
  const weight_t GOLDEN_RATIO = 0.38197;

  weight_t f_left = INF;
  weight_t f_mid = ComputeFunction(x + d * mid);
  weight_t f_right = INF;

  cerr << "Golden section set initial bracket [" << left << ", " 
       << mid << "(" << f_mid << "), " << right << "]" << endl;

  int iteration = 1;
  while (iteration < MAX_ITERATIONS && (right - left)/max(fabs(right),fabs(left)) > FINAL_BRACKET_SIZE) { 
    weight_t trial = (right - mid > mid - left) ?
      right - GOLDEN_RATIO * (right - left) :
      left + GOLDEN_RATIO * (right - left);

    vector<weight_t> testPoint = x + d * trial;
    weight_t f_trial = ComputeFunction (testPoint);

    if (f_trial <= f_mid) {
      if (trial >= mid) { 
	left = mid; 
	f_left = f_mid;
      }
      else {
	right = mid;
	f_right = f_mid;
      }
      mid = trial; 
      f_mid = f_trial;
    } 
    else {
      if (trial < mid) { 
	left = trial; 
	f_left = f_trial;
      }
      else {
	right = trial;
	f_right = f_trial;
      }
    }

    cerr << "Golden section iteration " << iteration << ": [" << setprecision(6) << left;
    if (f_left != INF)
      cerr << "(" << f_left << ")";
    cerr << ", " << mid << "(" << f_mid << "), " << right;
    if (f_right != INF)
      cerr << "(" << f_right << ")";
    cerr << "]" << endl;

    iteration++;
  }

  return mid;
}


static inline weight_t maxw(weight_t a, weight_t b) {
  return ((a > b) ? a : b);
}

/* Normally etaPlus is 1.2, etaMinus is 0.5 */
void Minimizer::rProp(vector<weight_t>& x, vector<weight_t>& learningRate,
		      weight_t etaPlus, weight_t etaMinus, bool surge) {

  RpropState state;
  RpropState prevState;
  bool prevBacktrack = false;  //keeps track of whether we backtracked on the previous iteration
  int iteration = 1;

  //set up initial state
  state.holds.resize(x.size(), false);
  state.learningRates = learningRate;
  state.metaLearningRates.resize(x.size(), etaPlus);

  cerr << "Beginning Rprop optimization, etaPlus = " << etaPlus << ", etaMinus = " << etaMinus << endl;
  if (surge)
    cerr << "Rprop surge enabled" << endl;

  while (true) {
    vector<weight_t> g;
    weight_t functionValue = ComputeFunctionAndGradient(g, x);

    if (surge && iteration > 4 && functionValue > state.functionValue && !prevBacktrack) {
      //we were too aggressive and the function value increased
      //go back to the previous point and reset all meta-learning rates
      cerr << "Rprop: function value increased to " << functionValue 
	   << ", backtracking and resetting meta-learning rates" << endl;
      x = state.x;
      state.holds = prevState.holds;
      state.learningRates = prevState.learningRates;
      state.metaLearningRates.clear();
      state.metaLearningRates.resize(x.size(), etaPlus);
      prevBacktrack = true;
    }
    else {
      prevBacktrack = false;
      prevState = state;
      state.x = x;
      state.functionValue = functionValue;
      state.gradient = g;
    }
    
    /* update the learning rates */
    if (iteration == 1) prevState.gradient = state.gradient;
    int advancers = 0;
    int decliners = 0;
    int holds = 0;
    for (int i = 0; i < x.size(); i++) {
      if (state.holds[i]) {
	//don't update learning rate, but release hold for next iteration
	state.holds[i] = false;
	holds++;
      }
      else if ( sign(prevState.gradient[i]) == sign(state.gradient[i]) ) {
	state.learningRates[i] *= state.metaLearningRates[i];
	advancers++;
	if (surge)
	  state.metaLearningRates[i] *= etaPlus;  //increase meta-learning rate
      }
      else {
	state.learningRates[i] *= etaMinus;
	decliners++;
	state.holds[i] = true;  //do not update learning rate for i on the next iteration
	if (surge)
	  state.metaLearningRates[i] = etaPlus;  //reset meta-learning rate
      }
    }

    // update parameters
    for (int i=0; i<x.size(); i++) {
      x[i] = x[i] - sign(state.gradient[i]) * state.learningRates[i];
    }

    cerr << "Rprop iteration " << iteration
	 << ": fval = " << state.functionValue 
	 << ", ||grad|| = " << Norm(state.gradient) 
	 << ", ||x|| = " << Norm(x)
	 << ", ||lr|| = " << Norm(state.learningRates);
    if (surge)
      cerr << ", ||mlr|| = " << Norm (state.metaLearningRates);
    cerr << ", " << advancers << " inc, "
	 << decliners << " dec"
	 << endl;

    if (! Report(x, iteration, state.functionValue, Norm(state.learningRates)))
      break;

    iteration++;
  }
}


weight_t sign(weight_t x) {
  if (x >= 0) return 1;
  else return -1;
}

/* Standard linear algebra */

weight_t DotProduct (const vector<weight_t> &x, const vector<weight_t> &y){
  assert (x.size() == y.size());
  weight_t ret = 0;
  for (int i = 0; i < (int) x.size(); i++) ret += x[i] * y[i];
  return ret;
}

weight_t Norm (const vector<weight_t> &x){
  return sqrt(DotProduct (x,x));
}

vector<weight_t>& Scale(vector<weight_t>& x, const vector<weight_t> &scale) {
  assert (x.size() == scale.size());
  for (int i=0; i<x.size(); i++)
    x[i] /= scale[i];
  return x;
}

vector<weight_t>& Unscale(vector<weight_t>& x, const vector<weight_t> &scale) {
  assert (x.size() == scale.size());
  for (int i=0; i<x.size(); i++)
    x[i] *= scale[i];
  return x;
}

vector<weight_t> operator+(const vector<weight_t> &x, const vector<weight_t> &y){
  vector<weight_t> ret (x);
  for (int i = 0; i < (int) ret.size(); i++) ret[i] += y[i];
  return ret;
}

vector<weight_t> operator-(const vector<weight_t> &x, const vector<weight_t> &y){
  vector<weight_t> ret (x);
  for (int i = 0; i < (int) ret.size(); i++) ret[i] -= y[i];
  return ret;
}

vector<weight_t> &operator+= (vector<weight_t> &x, weight_t c){
  for (int i = 0; i < (int) x.size(); i++) x[i] += c;
  return x;
}

vector<weight_t> operator*(const vector<weight_t> &x, weight_t c){
  vector<weight_t> ret (x);
  for (int i = 0; i < (int) ret.size(); i++) ret[i] *= c;
  return ret;
}
