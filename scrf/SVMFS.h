/* A set of features based on a support vector machine */

#ifndef SVMFS_H
#define SVMFS_H

#include "TransitionFeature.h"
#include "NumericFS.h"
#include "svm.h"
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef MULTI
#include <mpi.h>
#endif

enum kernel_t {POLYNOMIAL, GAUSSIAN};
enum svmfs_t {POSITIONAL, CODON};
typedef struct svm_model svm_t;

class SVMFS : public NumericFS {
 public:
  svmfs_t svmfsType;
  kernel_t kernelType;
  int kernelOrder;
  int numBins;
  int length;  //length of window SVM looks at
  int offset;  //distance to left side of window from feature position
  weight_t sampleRate;
  vector<string> seqNames;
  vector<int> seqIDs;
  vector<weight_t> binLimits;
  svm_t* svm;
  vector<vector<int> > trainingTrue;   //a vector of sparse binary vectors
  vector<vector<int> > trainingDecoys;
  vector<vector<int> > testingTrue;
  vector<vector<int> > testingDecoys;
  vector<vector<int> > svLookup;
  vector<float> weights;  //full weight vector for w^T F classification
  vector<int> codonFeatures;  //maps trimer indices to codon feature number
  int numCodonFeatures;

  char* bins[3];
  float* coeffs[3];

  //for CV training
  weight_t cvAccuracy;
  bool cvComplete;
  weight_t C;

  SVMFS(int id, svmfs_t svmfsType, kernel_t kernelType, int kernelOrder, int numBins, int length,
	int offset, weight_t sampleRate, vector<string>& seqNames);
  virtual void setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq);

  weight_t calcDecisionValue(vector<int>& x);
  virtual bool buildInputVector(vector<int>& x, pos_t pos, strand_t strand);
 
  void setSVLookup();         //set lookup table for fast calcProb
  void setWeightVector();
  void setCodonFeatureTable();  //set table that maps DNA trimer indices to codon features
  weight_t cvSVM(); //train and get generalization accuracy
  weight_t getAccuracy(bool training);
  void setSVMTrainingParameters(struct svm_parameter& parameters, weight_t C, bool probabilities);
  void buildProblem(struct svm_problem& problem, bool training);
  void freeProblem(struct svm_problem& problem);
  svm_t* copySVM(svm_t* svm);

  void readSVM(ifstream& svmStream);
  void printSVM(ostream& svmStream);
  void makeInitialGuess();
  void addRandomDecoy(TransitionFeature& transition, bool training);
  void printAlignmentBlock(ostream& os, vector<int>& x); 

  //for debugging
  void testCalcDecisionValue();

#ifdef MULTI
  void sendSVM(int destinationID);
  void receiveSVM(int sourceID);
#endif
};

#endif
