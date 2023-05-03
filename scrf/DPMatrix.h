/* Holds the Viterbi, forward, and backward matrices, along with
   a matrix of links used by the dyanmic programming algorithms.
   Also keeps track of scaling factors needed for making the
   algorithms numerically stable. */

#ifndef DPMATRIX_H
#define DPMATRIX_H

#include "Globals.h"
#include "Segment.h"
#include <vector>

class DPMatrix {
 public:
  pos_t length;
  int numberOfStates;
  int numberOfImplicitStates;
  int numberOfTransitions;
  int numberOfTransitionTypes;

  weight_t** stateSFSSums;
  bool** allowedTransition;
  bool** allowedTransitionType;
  weight_t** transitionSFSSums;
  weight_t** alpha;
  weight_t** beta;
  weight_t** alphaStar;
  weight_t** betaStar;
  weight_t** alphaStarStar;
  weight_t** betaStarStar;
  
  weight_t** viterbi;
  weight_t** posteriors;
  weight_t** meaDecode;
  stateid_t** maxState;
  pos_t** maxStateDuration;
  vector<Segment> segments;  /* all the possible segments in the sequence */
  vector<Segment*>*** segmentStartLinks; /* links to segments by position for DP speed */
  vector<Segment*>*** segmentEndLinks;

  DPMatrix();
  void setSegments(vector<Segment>& newSegments);
  void setSegmentEndLinks();
  void setSegmentStartLinks();
  
  void freeAllMatrices();
  void freeEverything();
  void freeSegments();
  void freeSegmentEndLinkVectors();
  void freeSegmentStartLinkVectors();
  void freeSegmentEndLinks();
  void freeSegmentStartLinks();
};

#endif
