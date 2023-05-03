/////////////////////////////////////////////////////////////////
// Score.h
//
// Description:
//     Routines for numerical computations.
//
// Notes:
//     If EXACT is defined, make sure that it is also defined
//     when compiling Score.cc -- otherwise, there will be a
//     mismatch of types not detected at compile time!
/////////////////////////////////////////////////////////////////

#ifndef SCORE_H
#define SCORE_H

#include <cmath>

/////////////////////////////////////////////////////////////////
// Options
/////////////////////////////////////////////////////////////////

#ifdef SANITY
#define EXACT
#endif

#ifdef EXACT
typedef double FloatingPointType;
#else
typedef float FloatingPointType;
#endif

/////////////////////////////////////////////////////////////////
// class Score
/////////////////////////////////////////////////////////////////

class Score {

  static FloatingPointType powers2_base[71];
  static FloatingPointType *powers2;

  FloatingPointType mantissa;
  int exponent;

 public:

  static Score ZERO;
  static Score ONE;

  Score();
  Score (bool unused);
  Score (const Score &rhs);
  Score (FloatingPointType m, int e);
  static void Adjust (FloatingPointType &m, int &e);
  static const Score FromFloat (FloatingPointType x);
  FloatingPointType ToFloat() const;

  Score &operator= (const Score &rhs);
  bool operator== (const Score &rhs) const;
  bool operator< (const Score &rhs) const;

  const Score operator+ (const Score &rhs) const;
  const Score operator* (const Score &rhs) const;
  const Score operator/ (const Score &rhs) const;
  
  Score &operator+= (const Score &rhs);
  Score &operator*= (const Score &rhs);
  Score &operator/= (const Score &rhs);
};

/////////////////////////////////////////////////////////////////
// Score::Score()
//
// Description:
//     Default constructor.  "value" is left uninitialized.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     None.
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline Score::Score(){}

/////////////////////////////////////////////////////////////////
// Score::Adjust()
//
// Description:
//     "Fix" mantissa/exponent pair.
//
// Arguments:
//     m -- a floating point value
//     e -- an integer value
//
// Returns:
//     a new Score object
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline void Score::Adjust (FloatingPointType &m, int &e){
  int offset;
  if (m < 1e-15 || m > 1e15){
#ifdef EXACT
    m = frexp (m, &offset);
#else
    m = frexpf (m, &offset);
#endif
    e += offset;
  }
}

/////////////////////////////////////////////////////////////////
// Score::Score()
//
// Description:
//     Copy constructor.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     None.
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline Score::Score(const Score &rhs) : 
  mantissa (rhs.mantissa), exponent (rhs.exponent) {}

/////////////////////////////////////////////////////////////////
// Score::Score()
//
// Description:
//     Create object using a particular internal representation,
//     after adjustment.
//
// Arguments:
//     m -- a floating point value
//     e -- an integer value
//
// Returns:
//     None.
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline Score::Score (FloatingPointType m, int e) : mantissa (m), exponent (e) {
  Adjust (mantissa, exponent);
}

/////////////////////////////////////////////////////////////////
// Score::FromFloat()
//
// Description:
//     Convert from FloatingPointType to Score object.
//
// Arguments:
//     x -- a floating point value
//
// Returns:
//     a new Score object
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline const Score Score::FromFloat (FloatingPointType x){
  return Score (x, 0);
}

/////////////////////////////////////////////////////////////////
// Score::ToFloat()
//
// Description:
//     Convert from Score object to FloatingPointType.
//
// Arguments:
//     None.
//
// Returns:
//     a FloatingPointType value
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline FloatingPointType Score::ToFloat() const {
#ifdef EXACT
  return ldexp(mantissa, exponent);
#else
  return ldexpf(mantissa, exponent);
#endif
}


/////////////////////////////////////////////////////////////////
// Score::operator=()
//
// Description:
//     Set current object to be equal to "rhs".
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     current object
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline Score &Score::operator= (const Score &rhs){
  if (this != &rhs){
    mantissa = rhs.mantissa;
    exponent = rhs.exponent;
  }

  return *this;
}

/////////////////////////////////////////////////////////////////
// Score::operator==()
//
// Description:
//     Test two Score values for equality.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     true or false
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline bool Score::operator== (const Score &rhs) const {
  return (mantissa == rhs.mantissa && exponent == rhs.exponent);
}

/////////////////////////////////////////////////////////////////
// Score::operator+()
//
// Description:
//     Add two Score values.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     sum of two Score values
//
// Environment:
//     powers2 matrix must have been previous initialized.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline const Score Score::operator+ (const Score &rhs) const {
  return Score (mantissa + rhs.mantissa * powers2[rhs.exponent - exponent], exponent);
}

/////////////////////////////////////////////////////////////////
// Score::operator*()
//
// Description:
//     Multiply two Score values.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     product of two Score values
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline const Score Score::operator* (const Score &rhs) const {
  return Score (mantissa * rhs.mantissa, exponent + rhs.exponent);
}

/////////////////////////////////////////////////////////////////
// Score::operator/()
//
// Description:
//     Divide two Score values.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     quotient of two Score values
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline const Score Score::operator/ (const Score &rhs) const {
  return Score (mantissa / rhs.mantissa, exponent - rhs.exponent);
}

/////////////////////////////////////////////////////////////////
// Score::operator+=()
//
// Description:
//     Add two Score values and stores in current object.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     current object
//
// Environment:
//     powers2 matrix must have been previous initialized.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline Score &Score::operator+= (const Score &rhs){
  mantissa += rhs.mantissa * powers2[rhs.exponent - exponent];
  Adjust (mantissa, exponent);
  return *this;
}

/////////////////////////////////////////////////////////////////
// Score::operator*=()
//
// Description:
//     Multiply two Score values and stores in current object.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     current object
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline Score &Score::operator*= (const Score &rhs){
  mantissa *= rhs.mantissa;
  exponent += rhs.exponent;
  Adjust (mantissa, exponent);
  return *this;
}

/////////////////////////////////////////////////////////////////
// Score::operator/=()
//
// Description:
//     Divide two Score values and stores in current object.
//
// Arguments:
//     rhs -- a Score
//
// Returns:
//     current object
//
// Environment:
//     None.
//
// Side-effects:
//     None.
/////////////////////////////////////////////////////////////////

inline Score &Score::operator/= (const Score &rhs){
  mantissa /= rhs.mantissa;
  exponent -= rhs.exponent;
  Adjust (mantissa, exponent);
  return *this;
}

#endif
