/////////////////////////////////////////////////////////////////
// Score.cc
//
// Description:
//     Routines for numerical computations.
/////////////////////////////////////////////////////////////////

#include "Score.h"

Score unknown (true);
Score Score::ZERO = Score::FromFloat(0);
Score Score::ONE = Score::FromFloat(1);
FloatingPointType Score::powers2_base[71];
FloatingPointType *Score::powers2;

/////////////////////////////////////////////////////////////////
// Score::Score()
//
// Description:
//     Static constructor to build powers2[] table.
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
//     Initializes the powers2 matrix.
/////////////////////////////////////////////////////////////////

Score::Score (bool unused){
  powers2 = powers2_base + 35;
  for (int i = -35; i <= 35; i++){
#ifdef EXACT
    powers2[i] = exp2(i);
#else
    powers2[i] = exp2f(i);
#endif
  }
}

