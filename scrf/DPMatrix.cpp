#include "DPMatrix.h"

DPMatrix::DPMatrix() {
  stateSFSSums = NULL;
  allowedTransition = NULL;
  allowedTransitionType = NULL;
  transitionSFSSums = NULL;
  alpha = NULL;
  beta = NULL; 
  alphaStar = NULL;
  betaStar = NULL;
  alphaStarStar = NULL;
  betaStarStar = NULL;

  viterbi = NULL;
  posteriors = NULL;
  meaDecode = NULL;
  maxState = NULL;
  maxStateDuration = NULL;

  segmentStartLinks = NULL;
  segmentEndLinks = NULL;
}

void DPMatrix::freeSegmentEndLinks() {
  if (segmentEndLinks != NULL) {
    freeSegmentEndLinkVectors();
    for (int i=0; i<numberOfStates; i++)
      delete[] segmentEndLinks[i];
    delete[] segmentEndLinks;
    segmentEndLinks = NULL;
  }
}


void DPMatrix::freeSegmentStartLinks() {
  if (segmentStartLinks != NULL) {
    freeSegmentStartLinkVectors();
    for (int i=0; i<numberOfStates; i++)
      delete[] segmentStartLinks[i];
    delete[] segmentStartLinks;
    segmentStartLinks = NULL;
  }  
}

void DPMatrix::freeAllMatrices() {
 #pragma omp sections
  {
  #pragma omp section
  if (stateSFSSums != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] stateSFSSums[j];
    delete[] stateSFSSums;
    stateSFSSums = NULL;
  }
  #pragma omp section
  if (allowedTransition != NULL) {
    delete[] allowedTransition;
    allowedTransition = NULL;
  }
  #pragma omp section
  if (allowedTransitionType != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] allowedTransitionType[j];
    delete[] allowedTransitionType;
    allowedTransitionType = NULL;
  }
  #pragma omp section
  if (transitionSFSSums != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] transitionSFSSums[j];
    delete[] transitionSFSSums;
    transitionSFSSums = NULL;
  }
  #pragma omp section
  if (alpha != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] alpha[j];
    delete[] alpha;
    alpha = NULL;
  }
  #pragma omp section
  if (beta != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] beta[j];
    delete[] beta;
    beta = NULL;
  }
  #pragma omp section
  if (alphaStar != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] alphaStar[j];
    delete[] alphaStar;
    alphaStar = NULL;
  }
  #pragma omp section
  if (betaStar != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] betaStar[j];
    delete[] betaStar;
    betaStar = NULL;
  }
  #pragma omp section
  if (alphaStarStar != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] alphaStarStar[j];
    delete[] alphaStarStar;
    alphaStarStar = NULL;
  }
  #pragma omp section
  if (betaStarStar != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] betaStarStar[j];
    delete[] betaStarStar;
    betaStarStar = NULL;
  }
  #pragma omp section
  if (viterbi != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] viterbi[j];
    delete[] viterbi;
    viterbi = NULL;
  }
  #pragma omp section
  if (posteriors != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] posteriors[j];
    delete[] posteriors;
    posteriors = NULL;
  }
  #pragma omp section
  if (meaDecode != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] meaDecode[j];
    delete[] meaDecode;
    meaDecode = NULL;
  }
  #pragma omp section
  if (maxState != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] maxState[j];
    delete[] maxState;
    maxState = NULL;
  }
  #pragma omp section
  if (maxStateDuration != NULL) {
    for (pos_t j=0; j<length; j++)
      delete[] maxStateDuration[j];
    delete[] maxStateDuration;
    maxStateDuration = NULL;
  }
  #pragma omp section 
  {
    freeSegmentEndLinks();
  }
  #pragma omp section
  {
    freeSegmentStartLinks();
  }
 }
}

void DPMatrix::freeEverything() {
  time_t startTime = time(NULL);

  freeAllMatrices();
  freeSegments();

  //  cerr << "freeEverything took " << (time(NULL) - startTime) << " seconds" << endl;
}

void DPMatrix::freeSegments() {
  freeSegmentEndLinkVectors();
  freeSegmentStartLinkVectors();

  /* swap the segments vector with an empty vector to really free the memory */
  vector<Segment>().swap(segments);
}


void DPMatrix::freeSegmentEndLinkVectors() {
  if (segmentEndLinks == NULL) return;

  for (int i=0; i<segments.size(); i++) {
    if (segmentEndLinks[segments[i].state][segments[i].end] != NULL) {
      delete segmentEndLinks[segments[i].state][segments[i].end];
      segmentEndLinks[segments[i].state][segments[i].end] = NULL;
    }
  }
}

void DPMatrix::freeSegmentStartLinkVectors() {
  if (segmentStartLinks == NULL) return;

  for (int i=0; i<segments.size(); i++) {
    if (segmentStartLinks[segments[i].state][segments[i].start] != NULL) {
      delete segmentStartLinks[segments[i].state][segments[i].start];
      segmentStartLinks[segments[i].state][segments[i].start] = NULL;
    }
  }
}

void DPMatrix::setSegments(vector<Segment>& newSegments) {
  /* clean up links from previous segment list */
  freeSegmentEndLinkVectors();
  freeSegmentStartLinkVectors();

  /* copy over new list of segments */
  segments = newSegments;
}

void DPMatrix::setSegmentStartLinks() {
  if (numberOfImplicitStates == numberOfStates)
    return;  //no segments for regular CRF model

  if (segmentStartLinks == NULL) {
    /* allocate matrix of links */
    segmentStartLinks = new vector<Segment*>**[numberOfStates];
    for (int i=0; i<numberOfStates; i++) {
      segmentStartLinks[i] = new vector<Segment*>*[length];
      for (pos_t j=0; j<length; j++)
	segmentStartLinks[i][j] = NULL; /* initialize */
    }
  }
  else {
    freeSegmentStartLinkVectors();
  }
  
  for (int i=0; i<segments.size(); i++) {
    stateid_t stateid = segments[i].state;
    pos_t start = segments[i].start;
    if (segmentStartLinks[stateid][start] == NULL) 
      segmentStartLinks[stateid][start] = new vector<Segment*>();
    segmentStartLinks[stateid][start]->push_back( &(segments[i]) );
  }
}

void DPMatrix::setSegmentEndLinks() {
  if (numberOfImplicitStates == numberOfStates)
    return;  //no segments for regular CRF model

  if (segmentEndLinks == NULL) {
    /* allocate matrix of links */
    segmentEndLinks = new vector<Segment*>**[numberOfStates];
    for (int i=0; i<numberOfStates; i++) {
      segmentEndLinks[i] = new vector<Segment*>*[length];
      for (pos_t j=0; j<length; j++)
	segmentEndLinks[i][j] = NULL; /* initialize */
    }
  }
  else {
    freeSegmentEndLinkVectors();
  }
  
  for (int i=0; i<segments.size(); i++) {
    stateid_t stateid = segments[i].state;
    pos_t end = segments[i].end;
    if (segmentEndLinks[stateid][end] == NULL) 
      segmentEndLinks[stateid][end] = new vector<Segment*>();
    segmentEndLinks[stateid][end]->push_back( &(segments[i]) );
  }
}
