/* The semi-Markov conditional random field */

#ifndef SEMICRF_H
#define SEMICRF_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include "State.h"
#include "StateType.h"
#include "TransitionFeature.h"
#include "TransitionType.h"
#include "DPMatrix.h"
#include "AlignmentSequence.h"
#include "NumericSequence.h"
#include "ESTSequence.h"
#include "Globals.h"
#include "MaskingFeature.h"
#include "DNAKmerFS.h"
#include "DNAKmerPairFS.h"
#include "DNAKmerArrayFS.h"
#include "DNAKmerPairArrayFS.h"
#include "SVMFS.h"
#include "ESTPositionFS.h"
#include "ESTTransitionFS.h"
#include "Segmentation.h"
#include "Inline.h"

using namespace std;

class Segmentation;

class SemiCRF {
 public:
  vector<State> states;        /* holds the states of the sCRF.  This vector is
				  sorted so that explicit-length states all come
				  before any implicit-length state. */
  int numberOfImplicitStates;
  int numberOfImplicitTypes;   /* number of state types describing at least
				  one implicit state */
  vector<StateType> stateTypes;
  vector<TransitionFeature> transitions;
  vector<weight_t> transitionWeights;
  vector<TransitionType> transitionTypes;
  DPMatrix dpMatrix;
  AlignmentSequence* alignment;
  NumericSequence* numericSeq;
  ESTSequence* estSeq;
  vector<SequenceFeatureSet*> sequenceFeatureSets;
  vector<weight_t*> parameterMap;   /* keeps track of locations of all model parameters */

  bool** allowedLookup;  //a table that says what transitions are allowed at given k-mer
  int allowedKmerLength;
  int allowedKmerOffset;

  //constructors
  SemiCRF(string parameterFile);
  SemiCRF(const SemiCRF& other);
  ~SemiCRF();

  //some utility functions
  stateid_t getStateID(string s);
  stateid_t getStateTypeID(string s);
  int getTransitionTypeID(string s);
  SVMFS* getSVMFSByID(int id);

  //for reading parameters
  void readParameters(string filename);
  void readStates(ifstream& parameterStream);
  void readTransitions(ifstream& parameterStream);
  void readLengths(ifstream& parameterStream);
  void readMaskingFeatures(ifstream& parameterStream);
  void readDNAKmerFeatures(ifstream& parameterStream);
  void readDNAKmerPairFeatures(ifstream& parameterStream);
  void readDNAKmerArrayFeatures(ifstream& parameterStream);
  void readDNAKmerPairArrayFeatures(ifstream& parameterStream);
  void readSVMFeatures(ifstream& parameterStream);
  void readESTPositionFeatures(ifstream& parameterStream);
  void readESTTransitionFeatures(ifstream& parameterStream);

  //for writing parameters
  void writeParameters(string filename);
  void printParameters(ostream& parameterStream);
  void printStates(ostream& parameterStream);
  void printTransitions(ostream& parameterStream);
  void printLengths(ostream& parameterStream);
  void printMaskingFeatures(ostream& parameterStream);
  void printDNAKmerFeatures(ostream& parameterStream);
  void printDNAKmerPairFeatures(ostream& parameterStream);
  void printDNAKmerArrayFeatures(ostream& parameterStream);
  void printDNAKmerPairArrayFeatures(ostream& parameterStream);
  void printSVMFeatures(ostream& parameterStream);
  void printESTPositionFeatures(ostream& parameterStream);
  void printESTTransitionFeatures(ostream& parameterStream);

  //for dynamic programming
  void setAllShortcuts(); /* set all shortcuts for fast memory access */
  void precomputeForDP(Segmentation* fixedPath);
  void setAllowedLookup();
  void setAllowedTransitions();
  weight_t forward();  /* sets alpha matrix */ 
  weight_t backward(); /* sets beta matrix */
  weight_t viterbi(Segmentation& traceback);
  weight_t computePartitionFromAlpha();
  weight_t computePartitionFromBeta();
  void featureExpectations(vector<weight_t>& counts);
  weight_t computeAnnotationPosteriors(Segmentation& annotation, vector<weight_t>& posteriors, weight_t lambda);
  weight_t computeAlphaStarPPP(Segmentation& annotation, vector<weight_t>& posteriors);
  weight_t computeBetaStarPPP(Segmentation& annotation, vector<weight_t>& posteroirs);
  weight_t computeAlphaStar(Segmentation& annotation, weight_t lambda);
  weight_t computeBetaStar(Segmentation& annotation, weight_t lambda);
  weight_t computePosteriorSumNumeratorFromAlphaStar();
  weight_t computePosteriorSumNumeratorFromBetaStar(Segmentation& annotation, weight_t lambda);
  void posteriorSumNumeratorGradientPPP(vector<weight_t>& g, Segmentation& annotation, vector<weight_t>& posteriors);
  void posteriorSumNumeratorGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda);
  
  void computePositionalPosteriors();
  void computePositionalDifference(Segmentation& annotation, vector<weight_t>& posteriorDifference);
  weight_t meaDecode(Segmentation& traceback, weight_t kappa=1.0);
  weight_t mesaDecode(Segmentation& traceback, weight_t kappa=1.0);

  void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq);
  void setParameterMap();  /* builds the parameter map and sets the global IDs of all model parameters */
  void setParameters(const vector<weight_t>& newParameters);
  void getNeededKmersAndSequences(set<int>& dnaKs, set<int>& dnaSeqs, 
				  set<int>& dnaPairKs, set<pair<int,int> >& dnaPairSeqs, set<string>& seqNames);
  void getNeededSequences(set<string>& seqNames);
  void precomputeKmerIndices();

  NumericSequence* createNumericSequence(AlignmentSequence* alignment);

  void computeAlphaStarStar(Segmentation& annotation, Segmentation& correctAnnotation, 
			    vector<weight_t>& posteriorDifference, weight_t gamma, weight_t lambda);
  void computeBetaStarStar(Segmentation& annotation, Segmentation& correctAnnotation, 
			   vector<weight_t>& posteriorDifference, weight_t gamma, weight_t lambda);
  void computePosteriorDifference(Segmentation& annotation, vector<weight_t>& posteriorDifference);
  weight_t AEasyGradTerm(Segmentation& annotation, Segmentation& correctAnnotation,
			 vector<weight_t>& posteriorDifference, weight_t gamma, weight_t lambda);
  void AHardGradTerm(Segmentation& annotation, Segmentation& correctAnnotation, vector<weight_t>& posteriorDifference, 
		     weight_t gamma, weight_t lambda, vector<weight_t>& g);
  void getMaxPosteriorIncorrectLabels(Segmentation& annotation, Segmentation& MPIL);

  weight_t computeESS(Segmentation& annotation, weight_t lambda);
  weight_t computeESSAlphaStar(Segmentation& annotation, weight_t lambda);
  weight_t computeESSBetaStar(Segmentation& annotation, weight_t lambda);
  void essNumeratorGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda);
  weight_t computeESSNumeratorFromAlphaStar();
  weight_t computeESSNumeratorFromBetaStar(Segmentation& annotation, weight_t lambda);

  weight_t computeESA(Segmentation& annotation, weight_t lambda);
  weight_t computeESAAlphaStar(Segmentation& annotation, weight_t lambda);
  weight_t computeESABetaStar(Segmentation& annotation, weight_t lambda);
  void esaNumeratorGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda);
  void esaGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda, weight_t esa);
  weight_t computeESANumeratorFromAlphaStar();
  weight_t computeESANumeratorFromBetaStar(Segmentation& annotation, weight_t lambda);

  weight_t computeASA(Segmentation& annotation, weight_t lambda, weight_t gamma, weight_t& feMultiplier);
  weight_t computeASAAlphaStar(Segmentation& annotation, weight_t lambda, weight_t gamma);
  weight_t computeASABetaStar(Segmentation& annotation, weight_t lambda, weight_t gamma);
  void asaGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda, weight_t gamma, weight_t feMultiplier);
  weight_t computeASANumeratorFromAlphaStar();
  weight_t computeASANumeratorFromBetaStar(Segmentation& annotation, weight_t lambda, weight_t gamma);
};

#endif
