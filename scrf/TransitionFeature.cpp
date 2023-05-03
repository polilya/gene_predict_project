#include "TransitionFeature.h"

//make sure transition is allowed based on k-mer content
bool TransitionFeature::allowed(AlignmentSequence* alignment, pos_t pos) {
  bool allowed = true;

  if (toStart.size() == 0 && fromEnd.size() == 0)
    return true;
  else if (fromEnd.size() != 0 && toStart.size() == 0) {
    allowed = false;
    for (int i=0; i<fromEnd.size(); i++) {
      int kmerLength = strlen(fromEnd[i]);
      if (strncmp(alignment->sequenceArray[0] + pos - kmerLength, 
		  fromEnd[i], kmerLength) == 0)
	allowed = true;
    }
  }
  else if (toStart.size() != 0 && fromEnd.size() == 0) {
    allowed = false;
    for (int i=0; i<toStart.size(); i++) {
      int kmerLength = strlen(toStart[i]);
      if (strncmp(alignment->sequenceArray[0] + pos, 
		  toStart[i], kmerLength) == 0)
	allowed = true;
    }
  }
  else
    fatalError("Transitions with both fromEnd and toStart not supported\n");

  return allowed;
}
