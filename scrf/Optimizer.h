#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "SemiCRF.h"
#include "Minimizer.h"
#include "GTFEvaluation.h"
#include "WorkUnit.h"
#include <time.h>

enum crf_command_t {COMPUTE_TRAINING_OBJECTIVE, 
		    COMPUTE_TESTING_OBJECTIVE,
		    COMPUTE_OBJECTIVE_AND_GRADIENT,
		    PREDICT, 
		    SEND_TIMING_INFO, 
		    TRAINING_WORK_REASSIGNMENT, 
		    TESTING_WORK_REASSIGNMENT, 
		    SEND_LENGTHS,
		    OPTIMIZATION_COMPLETE};

enum objective_t {LIKELIHOOD, EA, A, PPP, ESS, SS, ESA, ASA};

enum decode_t {VITERBI, MEA, MESA};

enum parameter_group_t {TRANSITIONS, TRANSITION_SFS, STATE_SFS};
const int numParameterGroups = 3;

class Optimizer : public Minimizer {
 public:
  int id;        /* the ID of this process */
  int numProcs;  /* total number of processes involved in optimization */
  objective_t objective; /* determines objective function we will optimizime */
  decode_t decode; /* determines decoding algorithm */
  SemiCRF* scrf;

  vector<AlignmentSequence*> trainingAlignments;
  vector<AlignmentSequence*> testingAlignments;
  vector<NumericSequence*> trainingNumericSeqs;
  vector<NumericSequence*> testingNumericSeqs;
  vector<ESTSequence*> trainingESTSeqs;
  vector<ESTSequence*> testingESTSeqs;

  /* Segmentations are used for training and function evaluation */
  vector<Segmentation> trainingSegmentations;
  vector<Segmentation> testingSegmentations;

  /* GTFs are used for evaluating predictions */
  vector<GTF> trainingGTFs;
  vector<GTF> testingGTFs;

  vector<string> trainingAlignmentFiles;
  vector<string> trainingESTSeqFiles;
  vector<string> trainingGTFFiles;

  vector<string> testingAlignmentFiles;
  vector<string> testingESTSeqFiles;
  vector<string> testingGTFFiles;

  vector<WorkUnit> myTrainingWork;
  vector<WorkUnit> allTrainingWork;
  vector<WorkUnit> myTestingWork;
  vector<WorkUnit> allTestingWork;

  vector<weight_t> initialGuess;
  vector<weight_t> savedIG;

  weight_t trainingLength;  /* total length of all training alignments */
  weight_t testingLength;   /* total length of all testing alignments */
  weight_t trainingPositionsToOptimize; /* number of training positions annotated as a state we are optimizing */
  weight_t testingPositionsToOptimize;  /* number of testing positions annotated as a state we are optimizing */
  string gtfOutputDir; /* where to write predictions to disk while training */
  weight_t C; /* regularization coefficient */
  time_t startTime;
  time_t lastTime;
  bool scaling;  /* should we scale all features by the number of times they occur in the training data?
		    this is implemented in a sort of tricky way by dividing the parameters and/or
		    gradient values as appropriate */ 
  weight_t gamma;
  weight_t lambda;
  weight_t kappa;
  bool stochastic;
  int currentTrainingExample;  /* for stochastic optimization */

  vector<weight_t> scalingFactors;
  ofstream trainingLog;

  weight_t minHoldoutObjective;
  int minHoldoutIteration;
  vector<weight_t> minHoldoutParams;

  Optimizer(SemiCRF* scrf, string trainingList, string testingList, objective_t objective, decode_t decode);
  ~Optimizer();

  void ShuffleExamples();
  void optimize(bool newInitialGuess, weight_t pseudocount, string gtfOutputDir, bool scaling, weight_t gamma, 
		weight_t lambda, weight_t kappa, bool stochastic, weight_t C = 0);
  virtual bool Report (const vector<weight_t> &theta, int iteration, 
		       weight_t objective, weight_t step_length);
  virtual void ComputeGradient (vector<weight_t> &g, const vector<weight_t> &x);
  virtual weight_t ComputeFunction  (const vector<weight_t> &x);
  virtual weight_t ComputeFunctionAndGradient  (vector<weight_t> &g, const vector<weight_t> &x);

  void predict(const vector<weight_t>& x, GTFEvaluation& eval, GTFEvaluation& testingEval,
	       weight_t& trainingLabelAccuracy, weight_t& testingLabelAccuracy,
	       weight_t& trainingObjective, weight_t& testingObjective);
  void setWorkAssignment(vector<int>& assignment, bool training);
  void resetWorkAssignment();
  void printWorkStats(vector<int>& workSums, ostream& os);
  void makeInitialGuess(weight_t pseudocount, vector<weight_t>& initialGuess, bool newInitialGuess);
  void getLengths();
  weight_t computeObjective (const vector<weight_t> &x, bool training);

  weight_t computeLikelihood(int alignmentNumber, bool training);
  weight_t likelihoodAndGradient (int trainingExample, vector<weight_t> &g);

  weight_t computeEA(int alignmentNumber, bool training);
  weight_t posteriorSumAndGradient (int trainingExample, vector<weight_t> &g);

  weight_t computeA(int alignmentNumber, bool training);
  weight_t AAndGradient (int trainingExample, vector<weight_t>& g);

  weight_t computePPP(int alignmentNumber, bool training);
  weight_t PPPAndGradient (int trainingExample, vector<weight_t> &g);

  weight_t computeESS(int alignmentNumber, bool training);
  weight_t ESSAndGradient (int trainingExample, vector<weight_t> &g);

  weight_t computeESA(int alignmentNumber, bool training);
  weight_t ESAAndGradient (int trainingExample, vector<weight_t> &g);

  weight_t computeASA(int alignmentNumber, bool training);
  weight_t ASAAndGradient (int trainingExample, vector<weight_t> &g);

  weight_t regularizationTerm(const vector<weight_t> &x);
  void regularizationGradient(const vector<weight_t> &x, vector<weight_t> &g);
  void computeScalingFactors();
  void scaleGridSearch(vector<weight_t>& initialGuess);
  weight_t weightedLabelAccuracy(Segmentation& a, Segmentation& b);


#ifdef MULTI
  void sendCommand(crf_command_t command, int destinationID);
  void broadcastCommand(crf_command_t command);
#endif

  /* for debugging */
  void testGradient();
  void testForwardBackward();
  void testSVMScores();
};

#endif
