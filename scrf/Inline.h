#ifndef INLINE_H
#define INLINE_H

#include "Globals.h"
#include <vector>

/* Some global inline functions */

/* returns x^y for integers -- not efficient */
static inline int intpow(const int x, const int y) {
  int result = 1;
  for (int i=0; i<y; i++)
    result *= x;
  return result;
}

static inline weight_t logOfOnePlusExpX(const weight_t x) {
  if (x > (weight_t)500) {
    //cerr << "Precision lost, x = " << x << endl; 
    return x;
  }  
  return log((weight_t)1 + exp(x));
}

static inline weight_t logOfOneMinusExpX(const weight_t x) {
  if (x > (weight_t)500) {
    //cerr << "Precision lost, x = " << x << endl; 
    return x;
  }  
  return log((weight_t)1 - exp(x));
}

static inline weight_t logAdd (const weight_t x, const weight_t y) {
  if (x < LOG_ZERO_REDUCED)
    return y;
  if (y < LOG_ZERO_REDUCED)
    return x;

  if (x < y)
    return x + logOfOnePlusExpX(y-x);
  else
    return y + logOfOnePlusExpX(x-y);
}

static inline weight_t logSubtract(const weight_t x, const weight_t y) {
  if (x < y)
    fatalError("Cannot store negative number in log space");
  return x + logOfOneMinusExpX(y-x);
}

static inline weight_t logPlusEquals (weight_t& x, const weight_t y) {
  return x = logAdd(x, y);
}

static inline weight_t logMinusEquals (weight_t& x, const weight_t y) {
  return x = logSubtract(x, y);
}

static inline weight_t Q(weight_t x, weight_t gamma) {
  return 1.0 / (1.0 + exp(-1.0 * x * gamma));
}

static inline weight_t dQdx(weight_t x, weight_t gamma) {
  weight_t eToNegativeGammaX = exp(-1.0 * x * gamma);
  return gamma * eToNegativeGammaX / 
    ((1.0 + eToNegativeGammaX) * (1.0 + eToNegativeGammaX));
}

/*
static inline weight_t Q(weight_t x, weight_t gamma) {
  return x;
}

static inline weight_t dQdx(weight_t x, weight_t gamma) {
  return 1.0;
}
*/

/*
static inline weight_t Q(weight_t x, weight_t gamma) {
  if (x > 0.0) return 1.0;
  else return 0.0;
}

static inline weight_t dQdx(weight_t x, weight_t gamma) {
  return 0;
}
*/

#endif
