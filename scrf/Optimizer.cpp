#include "Optimizer.h"
#include "Globals.h"
#include "Scaler.h"
#include "RegularizationSelector.h"
#include <iomanip>
#include <set>
#include <sstream>

#ifdef MULTI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

Optimizer::~Optimizer() {
  for (int i=0; i<trainingAlignments.size(); i++)
    delete trainingAlignments[i];
  for (int i=0; i<testingAlignments.size(); i++)
    delete testingAlignments[i];
}

Optimizer::Optimizer(SemiCRF* scrf, string trainingList, string testingList,
		     objective_t objective, decode_t decode) {
  id = 0;
  numProcs = 1;
  this->scrf = scrf;
  this->objective = objective;
  this->decode = decode;

#ifdef MULTI
  MPI_Comm_rank(MPI_COMM_WORLD, &id);  /* get our id */
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs); /* get number of process */
#endif

  cerr << "Process " << id << " now beginning CRF optimization";
#ifdef _OPENMP
  cerr << " with " << omp_get_max_threads() << " thread(s)";
#endif
  cerr << endl;

  /* parse the training list */
  ifstream trainingListStream(trainingList.c_str());
  if (! trainingListStream.is_open())
    fatalError("Could not open training list");
  while (! trainingListStream.eof()) {
    string s;
    trainingListStream >> s;
    if (s.length() > 0)
      trainingAlignmentFiles.push_back(s);
    trainingListStream >> s;
    if (s.length() > 0)
      trainingESTSeqFiles.push_back(s);
    trainingListStream >> s;
    if (s.length() > 0)
      trainingGTFFiles.push_back(s);
  } 
  if (trainingAlignmentFiles.size() < numProcs)
    fatalError("Optimizer: More processes requested than training examples found");

  if (id == 0)
    cerr << "Assigning training work" << endl;

  /* Initially assign work naively, we'll (maybe) load balance after the first iteration */
  for (int i=0; i<trainingAlignmentFiles.size(); i++)
    allTrainingWork.push_back(WorkUnit(i, trainingAlignmentFiles[i], 
				       trainingGTFFiles[i], i % numProcs));
  vector<int> trainingAssignment;
  for (int i=id; i<trainingAlignmentFiles.size(); i += numProcs)
    trainingAssignment.push_back(i);
  setWorkAssignment(trainingAssignment, true);

  /* now parse the list of testing files */
  if (testingList != "") {
    ifstream testingListStream(testingList.c_str());
    if (! testingListStream.is_open())
      fatalError("Could not open testing list");
    while (! testingListStream.eof()) {
      string s;
      testingListStream >> s;
      if (s.length() > 0)
	testingAlignmentFiles.push_back(s);
      testingListStream >> s;
      if (s.length() > 0)
	testingESTSeqFiles.push_back(s);
      testingListStream >> s;
      if (s.length() > 0)
	testingGTFFiles.push_back(s);
    }
  }

  if (id == 0)
    cerr << "Assigning testing work" << endl;

  /* Initially assign work naively, we'll (maybe) load balance after the first iteration */
  for (int i=0; i<testingAlignmentFiles.size(); i++)
    allTestingWork.push_back(WorkUnit(i, testingAlignmentFiles[i], testingGTFFiles[i],
				      i % numProcs));

  vector<int> testingAssignment;
  for (int i=id; i<testingAlignmentFiles.size(); i += numProcs)
    testingAssignment.push_back(i);
  setWorkAssignment(testingAssignment, false);

  /* open training log for writing */
  trainingLog.open("train.log");
}


void Optimizer::optimize(bool newInitialGuess, weight_t pseudocount, string gtfOutputDir, bool scaling, 
			 weight_t gamma, weight_t lambda, weight_t kappa, bool stoch, weight_t C) {

  this->stochastic = stoch;
  this->gtfOutputDir = gtfOutputDir;
  this->scaling = scaling;
  this->gamma = gamma;
  this->lambda = lambda;
  this->kappa = kappa;
  this->C = C;
  currentTrainingExample = 0;

  if (id == 0) {
    cerr << "Optimizing, C = " << this->C << endl;
    if (objective == EA || objective == A || objective == ESS || objective == SS || objective == ESA || objective == ASA)
      cerr << "Lambda = " << lambda << endl;
    if (objective == A || objective == SS)
      cerr << "Gamma = " << gamma << endl;
    if (decode == MEA || decode == MESA)
      cerr << "Kappa = " << kappa << endl;
  }

  if (newInitialGuess)
    makeInitialGuess(pseudocount, initialGuess, newInitialGuess);
  else {
    initialGuess.clear();
    initialGuess.resize(scrf->parameterMap.size());
    for (int i=0; i<scrf->parameterMap.size(); i++) {
      initialGuess[i] = *(scrf->parameterMap[i]);
    }
  }

  if (id == 0) {
    scrf->writeParameters("initial_parameters.txt");
  }

  savedIG = initialGuess;

  if (scaling) {
    computeScalingFactors();
    
    //adjust initial guess
    if (id == 0) {
      for (int i=0; i<initialGuess.size(); i++)
	initialGuess[i] *= scalingFactors[i];
      scrf->setParameters(initialGuess);
    }
  }

  if (id == 0) {    
    getLengths();

    //tests for debugging
    //testSVMScores();
    //testForwardBackward();
    //testGradient();

    trainingLog << "Total training alignment lengths are " << trainingLength << endl;
    trainingLog << "Total testing alignment lengths are " << testingLength << endl;
    trainingLog << "Model has " << scrf->parameterMap.size() << " total parameters" << endl;
    trainingLog << "Training positions to optimize: " << trainingPositionsToOptimize << endl;
    trainingLog << "Testing positions to optimize: " << testingPositionsToOptimize << endl;

   
      // If applicable, do grid search for scaling factors to adjust SVM initial guess
      for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
	if (scrf->sequenceFeatureSets[i]->type == SVM) {
	  scaleGridSearch(initialGuess);
	  break;
	}
      }
    


    startTime = time(NULL);
    trainingLog << "Beginning optimization procedure: " << numProcs << " processes" << endl;
    trainingLog << endl;

    /* make initial predictions */
    GTFEvaluation trainingEval;
    GTFEvaluation testingEval;
    weight_t trainingLabelAccuracy;
    weight_t testingLabelAccuracy;
    weight_t trainingObjective = 0;
    weight_t testingObjective = 0;
    predict(initialGuess, trainingEval, testingEval, trainingLabelAccuracy, testingLabelAccuracy,
	    trainingObjective, testingObjective);
    trainingLabelAccuracy /= trainingLength;
    testingLabelAccuracy /= testingLength;

    trainingLog << "Objective on the training data is " 
		<< -trainingObjective << " (" << -(trainingObjective-regularizationTerm(initialGuess)) 
		<< " without regularization)" << endl;
    trainingLog << "Objective on the testing data is " << -testingObjective << " (" 
		<< -(testingObjective-regularizationTerm(initialGuess)) 
		<< " without regularization)" << endl;
    trainingLog << endl;

    trainingLog << "Initial label accuracy on the training data is " << trainingLabelAccuracy << endl;
    trainingLog << "Initial label accuracy on the testing data is " << testingLabelAccuracy << endl;
    trainingLog << endl;

    trainingLog << "INITIAL ACCURACY ON TRAINING DATA" << endl;
    trainingLog << trainingEval << endl;
    trainingLog << endl;
    trainingLog << "INITIAL ACCURACY ON TESTING DATA" << endl;
    trainingLog << testingEval << endl;
    trainingLog << endl;
    
    lastTime = time(NULL);

    RegularizationSelector rs(this, initialGuess, trainingLog);
    vector<weight_t> zero(1, 0.0);
    vector<weight_t> one(1, 1.0);
    rs.GoldenSection(zero, one, -35, -20, -10);
    scrf->setParameters(rs.bestParams);
    scrf->writeParameters("final_parameters.txt");

    /*
    if (stochastic) {
      vector<weight_t> learningRate(initialGuess.size(), 1e-6);
      StochasticDescent(initialGuess, 0.1, learningRate, 0.9, 5*trainingAlignments.size(), 1000);
    }
    else {
      vector<weight_t> learningRate(initialGuess.size());
      for (int i=0; i<initialGuess.size(); i++)
	learningRate[i] = max(1e-4, 1e-3 * fabs(initialGuess[i]));
      rProp(initialGuess, learningRate);
    
      //LBFGS(initialGuess);
    }
    */
#ifdef MULTI
    broadcastCommand(OPTIMIZATION_COMPLETE);
#endif
  }
#ifdef MULTI
  else { /* listen for commands */
    while (true) {
      crf_command_t command;      
      MPI_Bcast(&command, 1, MPI_INT, 0, MPI_COMM_WORLD);  /* receive command */

      if (command == COMPUTE_TRAINING_OBJECTIVE) {
	/* receive new parameter values to use for processing this command */
	vector<weight_t> newParameters(scrf->parameterMap.size());
	MPI_Bcast(&newParameters.front(), newParameters.size(), MPI_WEIGHT_T, 0, MPI_COMM_WORLD);
	computeObjective(newParameters, true); 
      }
      else if (command == COMPUTE_TESTING_OBJECTIVE) {
	/* receive new parameter values to use for processing this command */
	vector<weight_t> newParameters(scrf->parameterMap.size());
	MPI_Bcast(&newParameters.front(), newParameters.size(), MPI_WEIGHT_T, 0, MPI_COMM_WORLD);
	computeObjective(newParameters, false); 
      }
      else if (command == COMPUTE_OBJECTIVE_AND_GRADIENT) {
	/* receive new parameter values to use for processing this command */
	vector<weight_t> newParameters(scrf->parameterMap.size());
	MPI_Bcast(&newParameters.front(), newParameters.size(), MPI_WEIGHT_T, 0, MPI_COMM_WORLD);
	vector<weight_t> gradient;
	ComputeFunctionAndGradient(gradient, newParameters); 
      }
      else if (command == PREDICT) {
	/* receive new parameter values to use for processing this command */
	vector<weight_t> newParameters(scrf->parameterMap.size());
	MPI_Bcast(&newParameters.front(), newParameters.size(), MPI_WEIGHT_T, 0, MPI_COMM_WORLD);

	/* make predictions on our data */
	GTFEvaluation trainingEval;
	GTFEvaluation testingEval;
	weight_t trainingLabelAccuracy;
	weight_t testingLabelAccuracy;
	weight_t trainingObjective = 0;
	weight_t testingObjective = 0;
	predict(newParameters, trainingEval, testingEval, trainingLabelAccuracy, testingLabelAccuracy,
		trainingObjective, testingObjective);
	/* send predictions to root */
	vector<int> packedTrainingEval;
	vector<int> packedTestingEval;
	trainingEval.pack(packedTrainingEval);
	testingEval.pack(packedTestingEval);
	MPI_Reduce(&packedTrainingEval.front(), NULL, packedTrainingEval.size(), 
		   MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&packedTestingEval.front(), NULL, packedTestingEval.size(), 
		   MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&trainingLabelAccuracy, NULL, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&testingLabelAccuracy, NULL, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      else if (command == SEND_TIMING_INFO) {
	vector<int> myTimingInfo(allTrainingWork.size(), 0);
	for (int i=0; i<myTrainingWork.size(); i++) {
	  myTimingInfo[myTrainingWork[i].id] = myTrainingWork[i].elapsedTime;
	  myTrainingWork[i].elapsedTime = 0;
	}
	MPI_Reduce(&myTimingInfo.front(), NULL, myTimingInfo.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      else if (command == TRAINING_WORK_REASSIGNMENT) {
	vector<int> assignment(allTrainingWork.size(), -1);
	MPI_Recv(&assignment.front(), assignment.size(), MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
	setWorkAssignment(assignment, true);
      }
      else if (command == TESTING_WORK_REASSIGNMENT) {
	vector<int> assignment(allTestingWork.size(), -1);
	MPI_Recv(&assignment.front(), assignment.size(), MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
	setWorkAssignment(assignment, false);
      }
      else if (command == SEND_LENGTHS) {
	getLengths();
      }
      else if (command == OPTIMIZATION_COMPLETE) {
	return;
      }
      else {
	cerr << "Process " << id << " received unrecognized command type " << command << endl;
	fatalError("Unrecognized command during distributed training");
      }
    }
  }
  #endif
}


bool Optimizer::Report (const vector<weight_t> &theta, int iteration, 
			weight_t objective, weight_t step_length)  {
  // write current parameters to disk
  vector<weight_t> thetaPrime = theta;
  if (scaling) {
    for (int i=0; i<thetaPrime.size(); i++)
      thetaPrime[i] /= scalingFactors[i];
  }
  scrf->setParameters(thetaPrime);
  ostringstream oss;
  oss << "iteration" << iteration << ".parameters.txt";
  scrf->writeParameters(oss.str());


  //compute objective on the holdout set, not counting regularization
  weight_t holdoutObjective = computeObjective(theta, false) - regularizationTerm(theta);
  if (iteration == 1 || holdoutObjective < minHoldoutObjective) {
    minHoldoutObjective = holdoutObjective;
    minHoldoutParams = theta;
    minHoldoutIteration = iteration;
  }
  cerr << "Holdout objective is " << holdoutObjective << endl;

  time_t totalTime = time(NULL) - startTime; 
  time_t iterTime = time(NULL) - lastTime;
  lastTime = time(NULL);

  trainingLog << "Iteration " << iteration << " complete" << endl;
  trainingLog << "Training objective is " << objective << ", holdout objective is " << holdoutObjective << endl;
  if (stochastic) {
    trainingLog  << "Learning rate is " << step_length << endl;
  }
  trainingLog << "Time elapsed during this iteration: " << iterTime << " seconds" << endl;
  trainingLog << "Total time elapsed: " << totalTime << " seconds" << endl;
  trainingLog << endl; 
  trainingLog.flush(); 

  if (iteration - minHoldoutIteration >= 5) {  //looks like we're overfitting, optimization is complete 
    if (testingGTFs.size() != 0) {
      // predict and evaluate on both the training and testing data sets using best parameters
      GTFEvaluation trainingEval;
      GTFEvaluation testingEval;
      weight_t trainingLabelAccuracy;
      weight_t testingLabelAccuracy;
      weight_t trainingObjective = 0;
      weight_t testingObjective = 0;
      predict(minHoldoutParams, trainingEval, testingEval, trainingLabelAccuracy, testingLabelAccuracy,
	      trainingObjective, testingObjective);
      trainingLabelAccuracy /= trainingLength;
      testingLabelAccuracy /= testingLength;
      
      trainingLog << "Objective on the training data is " << -trainingObjective << " (" 
		  << -(trainingObjective-regularizationTerm(theta)) 
		  << " without regularization)" << endl;
      trainingLog << "Objective on the testing data is " << -testingObjective << " (" 
		  << -(testingObjective-regularizationTerm(theta)) 
		  << " without regularization)" << endl;
      trainingLog << endl;
      
      trainingLog << "Label accuracy on the training data is " << trainingLabelAccuracy << endl;
      trainingLog << "Label accuracy on the testing data is " << testingLabelAccuracy << endl << endl;
      
      trainingLog << "ACCURACY ON TRAINING DATA" << endl;
      trainingLog << trainingEval;
      trainingLog << endl;
      
      trainingLog << "ACCURACY ON TESTING DATA" << endl;
      trainingLog << testingEval;
      trainingLog << endl;
    }
    return false;
  }
  else
    return true;  //continue optimization
}

void Optimizer::predict(const vector<weight_t>& x, GTFEvaluation& trainingEval, GTFEvaluation& testingEval,
			weight_t& trainingLabelAccuracy, weight_t& testingLabelAccuracy,
			weight_t& trainingObjective, weight_t& testingObjective) {

  time_t startTime;
  if (id == 0)
    startTime = time(NULL);


  vector<weight_t> xPrime = x;

  if (id == 0 && scaling) {
    for (int i=0; i<x.size(); i++)
      xPrime[i] /= scalingFactors[i];
  }

#ifdef MULTI
  if (id == 0) {
    /* broadcast command to make predictions */
    crf_command_t command = PREDICT;
    MPI_Bcast(&command, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* broadcast parameters to be used */
    MPI_Bcast(&xPrime.front(), xPrime.size(), MPI_WEIGHT_T, 0, MPI_COMM_WORLD);
  }    
#endif

  scrf->setParameters(xPrime);

  /* run predictions on the training data */
  trainingLabelAccuracy = (weight_t)0;
  for (int i=0; i<trainingAlignments.size(); i++) {
    scrf->setSequences(trainingAlignments[i], trainingNumericSeqs[i], trainingESTSeqs[i]);
    scrf->precomputeForDP(NULL);
    Segmentation prediction;
    if (decode == VITERBI)
      scrf->viterbi(prediction);
    else if (decode == MEA || decode == MESA) {
#pragma omp parallel sections
      {
#pragma omp section
	scrf->forward();
#pragma omp section	
	scrf->backward();
      }
      if (decode == MEA) {
	scrf->computePositionalPosteriors();
	scrf->meaDecode(prediction, kappa);
      }
      else if (decode == MESA)
	scrf->mesaDecode(prediction, kappa);
    }
    prediction.setLabels();
    trainingLabelAccuracy += weightedLabelAccuracy(trainingSegmentations[i], prediction);
    GTF predictionGTF(prediction, *scrf);
    GTFEvaluation performance(trainingGTFs[i], predictionGTF);
    trainingEval += performance;
  }

  /* run predictions on the testing data */
  testingLabelAccuracy = (weight_t)0;
  for (int i=0; i<testingAlignments.size(); i++) {
    scrf->setSequences(testingAlignments[i], testingNumericSeqs[i], testingESTSeqs[i]);
    scrf->precomputeForDP(NULL);
    Segmentation prediction;
    if (decode == VITERBI)
      scrf->viterbi(prediction);
    else if (decode == MEA || decode == MESA) {
#pragma omp parallel sections
      {
#pragma omp section
	scrf->forward();
#pragma omp section	
	scrf->backward();
      }
      if (decode == MEA) {
	scrf->computePositionalPosteriors();
	scrf->meaDecode(prediction, kappa);
      }
      else if (decode == MESA)
	scrf->mesaDecode(prediction, kappa);
    }
    prediction.setLabels();
    testingLabelAccuracy += weightedLabelAccuracy(testingSegmentations[i], prediction);
    GTF predictionGTF(prediction, *scrf);
    GTFEvaluation performance(testingGTFs[i], predictionGTF);
    testingEval += performance;
    
    if (gtfOutputDir != "NULL") {
      /* write predictions to disk */
      int testingIndex = id + i * numProcs;
      string testingGTFFile = gtfOutputDir + "/" + testingGTFFiles[testingIndex];
      if (testingGTFFile.substr(testingGTFFile.length() - 3, 3) != "gtf")
	fatalError("Testing GTF filename does not end with \"gtf\"");
      testingGTFFile.replace(testingGTFFile.length() - 3, 3, "pred.gtf");
      predictionGTF.write(testingGTFFile);
    }
  }

#ifdef MULTI
  if (id == 0) {
    vector<int> packedTrainingEval;
    vector<int> packedTestingEval;
    trainingEval.pack(packedTrainingEval);
    testingEval.pack(packedTestingEval);
    vector<int> packedCombinedTrainingEval(packedTrainingEval.size(), 0);
    vector<int> packedCombinedTestingEval(packedTestingEval.size(), 0);
    weight_t combinedTrainingLabelAccuracy;
    weight_t combinedTestingLabelAccuracy;

    /* receive accuracy reports from other processes and combine them with ours */
    MPI_Reduce(&packedTrainingEval.front(), &packedCombinedTrainingEval.front(), 
	       packedTrainingEval.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&packedTestingEval.front(), &packedCombinedTestingEval.front(), 
	       packedTestingEval.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&trainingLabelAccuracy, &combinedTrainingLabelAccuracy, 1, 
	       MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&testingLabelAccuracy, &combinedTestingLabelAccuracy, 1, 
	       MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);

    trainingEval.unpack(packedCombinedTrainingEval);
    testingEval.unpack(packedCombinedTestingEval);
    trainingLabelAccuracy = combinedTrainingLabelAccuracy;
    testingLabelAccuracy = combinedTestingLabelAccuracy;
  }
#endif

  if (id == 0)
    cerr << "predict() completed in " << (time(NULL) - startTime) << " seconds" << endl;
}

weight_t Optimizer::weightedLabelAccuracy(Segmentation& a, Segmentation& b) {
  if (a.length != b.length)
    fatalError ("Cannot compute weightedLabelAccuracy for segmentations of different lengths");

  weight_t accuracy = 0;

  for (int j=0; j<a.length; j++) {
    if (a.labels[j] == b.labels[j]) {
      if (scrf->states[a.labels[j]].optimize)
	accuracy += lambda;
      else
	accuracy += 1.0;
    }
  }

  return accuracy;
}

void Optimizer::ComputeGradient(vector<weight_t> &g, const vector<weight_t> &x) {
  ComputeFunctionAndGradient(g, x);
}

void Optimizer::ShuffleExamples() {
  int numTimes = 3;

  int size = trainingAlignments.size();
  assert(trainingNumericSeqs.size() == size);
  assert(trainingSegmentations.size() == size);
  assert(trainingGTFs.size() == size);
  assert(myTrainingWork.size() == size);
  
  for (int i = 0; i < numTimes; i++) {
    for (int j = 0; j < size; j++) {
      int newIndex = j + rand() % (size - j);

      vector<AlignmentSequence *>::iterator a1 =
	trainingAlignments.begin() + j;
      vector<AlignmentSequence *>::iterator a2 = 
	trainingAlignments.begin() + newIndex;
      swap(a1, a2);

      vector<NumericSequence*>::iterator b1 =
	trainingNumericSeqs.begin() + j;
      vector<NumericSequence*>::iterator b2 =
	trainingNumericSeqs.begin() + newIndex;
      swap(b1, b2);
      
      vector<Segmentation>::iterator c1 =
	trainingSegmentations.begin() + j;
      vector<Segmentation>::iterator c2 =
	trainingSegmentations.begin() + newIndex;
      swap(c1, c2);

      vector<GTF>::iterator d1 =
	trainingGTFs.begin() + j;
      vector<GTF>::iterator d2 =
	trainingGTFs.begin() + newIndex;    
      swap(d1, d2);

      vector<WorkUnit>::iterator f1 =
	myTrainingWork.begin() + j;
      vector<WorkUnit>::iterator f2 =
	myTrainingWork.begin() + newIndex;
      swap(f1, f2);
    }
  }

}

weight_t Optimizer::ComputeFunctionAndGradient (vector<weight_t> &g, const vector<weight_t> &x) {
  time_t startTime;
  if (id == 0)
    startTime = time(NULL);
  
  vector<weight_t> xPrime = x;

  if (id == 0 && scaling) {
    for (int i=0; i<x.size(); i++)
      xPrime[i] /= scalingFactors[i];
  }
  
#ifdef MULTI
  if (id == 0) {

    /* broadcast command to begin gradient calculation */
    crf_command_t command = COMPUTE_OBJECTIVE_AND_GRADIENT;
    MPI_Bcast(&command, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    /* send parameters to be used */
    MPI_Bcast(&xPrime.front(), xPrime.size(), MPI_WEIGHT_T, 0, MPI_COMM_WORLD);
  }
#endif

  scrf->setParameters(xPrime);

  g.clear();
  g.resize(scrf->parameterMap.size(), 0);
  weight_t functionVal = 0;
 
  if (objective == LIKELIHOOD || objective == EA || objective == A || objective == PPP 
      || objective == ESS || objective == ESA || objective == ASA) {
    for (int i=0; i<trainingAlignments.size(); i++) {

      /* If stochastic, only go through this loop once, for the desired
       * example, and then break */
      if (stochastic)
	  i = currentTrainingExample;

      myTrainingWork[i].startTime = time(NULL);
      
      vector<weight_t> gradientContribution(g.size());
      if (objective == LIKELIHOOD)
	functionVal += likelihoodAndGradient(i, gradientContribution);
      else if (objective == EA)
	functionVal += posteriorSumAndGradient(i, gradientContribution);
      else if (objective == A)
	functionVal += AAndGradient(i, gradientContribution);
      else if (objective == PPP)
	functionVal += PPPAndGradient(i, gradientContribution);
      else if (objective == ESS)
	functionVal += ESSAndGradient(i, gradientContribution);
      else if (objective == ESA)
	functionVal += ESAAndGradient(i, gradientContribution);
      else if (objective == ASA)
	functionVal += ASAAndGradient(i, gradientContribution);

      for (int j=0; j<g.size(); j++) {
	g[j] += gradientContribution[j];
      }
      
      myTrainingWork[i].endTime = time(NULL);
      myTrainingWork[i].elapsedTime += myTrainingWork[i].endTime - myTrainingWork[i].startTime;

      if (stochastic) break;
    }
     
#ifdef MULTI
    if (id == 0) {
      if (stochastic) {
	/* receive length contributions from other processes */
	MPI_Reduce(&trainingAlignments[currentTrainingExample]->length, &trainingLength, 
		   1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      }

      /* receive function value contributions from other processes and add them to ours */
      weight_t combinedFunctionVal;
      MPI_Reduce(&functionVal, &combinedFunctionVal, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
      functionVal = combinedFunctionVal;

      /* receive gradient contributions from other processes and add them to ours*/
      vector<weight_t> combinedGradient(g.size(), 0);
      MPI_Reduce(&g.front(), &combinedGradient.front(), g.size(), MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
      for (int i=0; i<g.size(); i++)
	g[i] = combinedGradient[i];
    }
    else {
      if (stochastic) {
	/* send the length of our current training example to the root process */
	MPI_Reduce(&trainingAlignments[currentTrainingExample]->length, NULL, 1,
		   MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      }

      /* send our function value contribution to the root process */
      MPI_Reduce(&functionVal, NULL, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);

      /* send our gradient contribution to the root process */
      MPI_Reduce(&g.front(), NULL, g.size(), MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif  
    if (id == 0) {
      //normalize by length
      functionVal /= trainingLength;
      for (int i=0; i<g.size(); i++)
	g[i] /= trainingLength;
    }
  }    

  if (id == 0) {
    if (scaling) {
      for (int i=0; i<g.size(); i++)
	g[i] /= scalingFactors[i];
    }    
     
    // add regularization term
    functionVal += regularizationTerm(x);
    vector<weight_t> regGradient;
    regularizationGradient(x, regGradient);
    for (int i=0; i<g.size(); i++)
      g[i] += regGradient[i];
  }
  
  if (id == 0) {
    cerr << "ComputeFunctionAndGradient() completed in " << time(NULL) - startTime << " seconds" << endl;
  }

  /* update the current training example */
  if (stochastic) {
    currentTrainingExample++;
    if (currentTrainingExample >= trainingAlignments.size()) {
      currentTrainingExample = 0;
      ShuffleExamples();
    }
  }

  return functionVal;
}

weight_t Optimizer::likelihoodAndGradient (int trainingExample, vector<weight_t>& g) {
    vector<weight_t> trainingExpectations(scrf->parameterMap.size());
    vector<weight_t> allExpectations(scrf->parameterMap.size());
    weight_t trainingPartition, allPartition;
    int threadID;

    scrf->setSequences(trainingAlignments[trainingExample], trainingNumericSeqs[trainingExample], trainingESTSeqs[trainingExample]);
    scrf->precomputeForDP(&(trainingSegmentations[trainingExample]));

#pragma omp parallel sections
 {
  #pragma omp section
   trainingPartition = scrf->forward();
  #pragma omp section
   scrf->backward();
 }

 scrf->featureExpectations(trainingExpectations);
 scrf->dpMatrix.freeEverything();
 scrf->precomputeForDP(NULL);

#pragma omp parallel sections
 {
  #pragma omp section
     allPartition = scrf->forward();
  #pragma omp section
     scrf->backward();
 }

 scrf->featureExpectations(allExpectations);

 for (int i=0; i<g.size(); i++) {
   /* optimize negative log likelihood */
   g[i] = (weight_t)-1 * (trainingExpectations[i] - allExpectations[i]);
 }
 
 return (weight_t)-1 * (trainingPartition - allPartition);
}


weight_t Optimizer::posteriorSumAndGradient (int trainingExample, vector<weight_t>& g) {
  scrf->setSequences(trainingAlignments[trainingExample], trainingNumericSeqs[trainingExample], trainingESTSeqs[trainingExample]);
  scrf->precomputeForDP(NULL);

#pragma omp parallel sections 
 {
   #pragma omp section
   scrf->forward();
   
   #pragma omp section
   scrf->backward();
 }

 vector<weight_t> posteriors;
 weight_t posteriorSum = scrf->computeAnnotationPosteriors(trainingSegmentations[trainingExample], 
							   posteriors, lambda);

#pragma omp parallel sections
 {
   #pragma omp section
   scrf->computeAlphaStar(trainingSegmentations[trainingExample], lambda);

   #pragma omp section
   scrf->computeBetaStar(trainingSegmentations[trainingExample], lambda);
 }

  vector<weight_t> partitionGradient(g.size());
  vector<weight_t> posteriorSumNumeratorGradient(g.size()); 
  scrf->featureExpectations(partitionGradient);
  scrf->posteriorSumNumeratorGradient(posteriorSumNumeratorGradient, trainingSegmentations[trainingExample], lambda);

  for (int i=0; i<g.size(); i++)
    /* optimize negative expected accuracy */
    g[i] = (weight_t)-1 * (posteriorSumNumeratorGradient[i] - posteriorSum * partitionGradient[i]);

  return (weight_t)-1 * posteriorSum;
}

weight_t Optimizer::ESSAndGradient (int trainingExample, vector<weight_t>& g) {
  scrf->setSequences(trainingAlignments[trainingExample], trainingNumericSeqs[trainingExample], trainingESTSeqs[trainingExample]);
  scrf->precomputeForDP(NULL);

#pragma omp parallel sections 
 {
   #pragma omp section
   scrf->forward();
   
   #pragma omp section
   scrf->backward();
 }

 weight_t ess = scrf->computeESS(trainingSegmentations[trainingExample], lambda);

#pragma omp parallel sections
 {
   #pragma omp section
   scrf->computeESSAlphaStar(trainingSegmentations[trainingExample], lambda);

   #pragma omp section
   scrf->computeESSBetaStar(trainingSegmentations[trainingExample], lambda);
 }

  vector<weight_t> partitionGradient(g.size());
  vector<weight_t> essNumeratorGradient(g.size()); 
  scrf->featureExpectations(partitionGradient);
  scrf->essNumeratorGradient(essNumeratorGradient, trainingSegmentations[trainingExample], lambda);

  for (int i=0; i<g.size(); i++) {
    /* optimize negative expected accuracy */
    g[i] = (weight_t)-1 * (essNumeratorGradient[i] - ess * partitionGradient[i]);
  }

  return (weight_t)-1 * ess;
}

weight_t Optimizer::ESAAndGradient (int trainingExample, vector<weight_t>& g) {
  scrf->setSequences(trainingAlignments[trainingExample], trainingNumericSeqs[trainingExample], trainingESTSeqs[trainingExample]);
  scrf->precomputeForDP(NULL);

#pragma omp parallel sections
 {
   #pragma omp section
   scrf->computeESAAlphaStar(trainingSegmentations[trainingExample], lambda);

   #pragma omp section
   scrf->computeESABetaStar(trainingSegmentations[trainingExample], lambda);
 }

 weight_t esa = exp(scrf->computeESANumeratorFromAlphaStar() - scrf->computePartitionFromAlpha());

 scrf->esaGradient(g, trainingSegmentations[trainingExample], lambda, esa);

 return esa;
}

weight_t Optimizer::ASAAndGradient (int trainingExample, vector<weight_t>& g) {
  scrf->setSequences(trainingAlignments[trainingExample], trainingNumericSeqs[trainingExample], trainingESTSeqs[trainingExample]);
  scrf->precomputeForDP(NULL);

#pragma omp parallel sections
  {
#pragma omp section
    scrf->forward();

#pragma omp section
    scrf->backward();
  }

#pragma omp parallel sections
  {
#pragma omp section
    scrf->computeASAAlphaStar(trainingSegmentations[trainingExample], lambda, gamma);

#pragma omp section
    scrf->computeASABetaStar(trainingSegmentations[trainingExample], lambda, gamma);
  }

  weight_t feMultiplier;
  weight_t asa = scrf->computeASA(trainingSegmentations[trainingExample], lambda, gamma, feMultiplier);
  
  scrf->asaGradient(g, trainingSegmentations[trainingExample], lambda, gamma, feMultiplier);

  return asa;
}

weight_t Optimizer::AAndGradient (int trainingExample, vector<weight_t>& g) {
  startTime = time(NULL);

  scrf->setSequences(trainingAlignments[trainingExample], trainingNumericSeqs[trainingExample], trainingESTSeqs[trainingExample]);
  scrf->precomputeForDP(NULL);


#pragma omp parallel sections
 {
  #pragma omp section
    scrf->forward();
  #pragma omp section
    scrf->backward();
 }


  scrf->computePositionalPosteriors();

  vector<weight_t> posteriorDifference;
  scrf->computePosteriorDifference(trainingSegmentations[trainingExample], posteriorDifference);

  weight_t aa = 0;
  for (pos_t j=0; j<posteriorDifference.size(); j++) {
    if (scrf->states[trainingSegmentations[trainingExample].labels[j]].optimize)
      aa += lambda * Q(posteriorDifference[j], gamma);
    else
      aa += Q(posteriorDifference[j], gamma);
  }

  weight_t partition = scrf->computePartitionFromAlpha();
  vector<weight_t> partitionGradient(g.size());
  scrf->featureExpectations(partitionGradient);

#pragma omp parallel sections
 {
   #pragma omp section
     scrf->computeAlphaStarStar(trainingSegmentations[trainingExample], trainingSegmentations[trainingExample],
				posteriorDifference, gamma, lambda);
   #pragma omp section  
     scrf->computeBetaStarStar(trainingSegmentations[trainingExample], trainingSegmentations[trainingExample],
			       posteriorDifference, gamma, lambda);
 }

  weight_t easyTermCorrect = scrf->AEasyGradTerm(trainingSegmentations[trainingExample], 
						 trainingSegmentations[trainingExample],
						 posteriorDifference, gamma, lambda);
  vector<weight_t> hardTermCorrect;
  scrf->AHardGradTerm(trainingSegmentations[trainingExample], trainingSegmentations[trainingExample],
		      posteriorDifference, gamma, lambda, hardTermCorrect);

  freeWeightMatrix(scrf->dpMatrix.alphaStarStar, scrf->alignment->length, scrf->states.size());
  scrf->dpMatrix.alphaStarStar = NULL;
  freeWeightMatrix(scrf->dpMatrix.betaStarStar, scrf->alignment->length, scrf->states.size());
  scrf->dpMatrix.betaStarStar = NULL;

  Segmentation MPIL;
  scrf->getMaxPosteriorIncorrectLabels(trainingSegmentations[trainingExample], MPIL);
  
#pragma omp parallel sections
 {
  #pragma omp section
    scrf->computeAlphaStarStar(MPIL, trainingSegmentations[trainingExample], posteriorDifference, gamma, lambda);
  #pragma omp section
    scrf->computeBetaStarStar(MPIL, trainingSegmentations[trainingExample], posteriorDifference, gamma, lambda);
 }  

  weight_t easyTermMaxIncorrect = scrf->AEasyGradTerm(MPIL, trainingSegmentations[trainingExample],
						      posteriorDifference, gamma, lambda);
  vector<weight_t> hardTermMaxIncorrect;
  scrf->AHardGradTerm(MPIL, trainingSegmentations[trainingExample], 
		      posteriorDifference, gamma, lambda, hardTermMaxIncorrect);

  for (int i=0; i<g.size(); i++)
    /* optimize negative expected accuracy */
    g[i] = (weight_t)-1 * (hardTermCorrect[i] - partitionGradient[i] * easyTermCorrect 
			   - hardTermMaxIncorrect[i] + partitionGradient[i] * easyTermMaxIncorrect);

  return (weight_t)-1 * aa;
}

weight_t Optimizer::PPPAndGradient (int trainingExample, vector<weight_t>& g) {
  scrf->setSequences(trainingAlignments[trainingExample], trainingNumericSeqs[trainingExample], trainingESTSeqs[trainingExample]);
  scrf->precomputeForDP(NULL);

#pragma omp parallel sections 
 {
   #pragma omp section
   scrf->forward();
   
   #pragma omp section
   scrf->backward();
 }

 vector<weight_t> posteriors;
 scrf->computeAnnotationPosteriors(trainingSegmentations[trainingExample], posteriors, lambda);
 weight_t ppp = 0;
 for (int i=0; i<posteriors.size(); i++)
   ppp += posteriors[i];

#pragma omp parallel sections
 {
   #pragma omp section
   scrf->computeAlphaStarPPP(trainingSegmentations[trainingExample], posteriors);

   #pragma omp section
   scrf->computeBetaStarPPP(trainingSegmentations[trainingExample], posteriors);
 }

  vector<weight_t> partitionGradient(g.size());
  vector<weight_t> posteriorSumNumeratorGradient(g.size()); 
  scrf->featureExpectations(partitionGradient);
  scrf->posteriorSumNumeratorGradientPPP(posteriorSumNumeratorGradient, trainingSegmentations[trainingExample], posteriors);

  for (int i=0; i<g.size(); i++)
    /* optimize negative sum of log posteriors */
    g[i] = (weight_t)-1 * (posteriorSumNumeratorGradient[i] - scrf->alignment->length * partitionGradient[i]);

  return (weight_t)-1 * ppp;
}

weight_t Optimizer::ComputeFunction (const vector<weight_t> &x) {
  return computeObjective(x, true);  //compute the objective on the training data
}

weight_t Optimizer::computeObjective (const vector<weight_t> &x, bool training) {

  time_t startTime;
  if (id == 0)
    startTime = time(NULL);

  vector<weight_t> xPrime = x;
  
  if (id == 0 && scaling) {
    for (int i=0; i<x.size(); i++)
      xPrime[i] /= scalingFactors[i];
  }
    
#ifdef MULTI
  if (id == 0) { 
    /* broadcast command to begin function evaluation */
    crf_command_t command;
    if (training)
      command = COMPUTE_TRAINING_OBJECTIVE;
    else
      command = COMPUTE_TESTING_OBJECTIVE;
    MPI_Bcast(&command, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    /* send parameters to be used */
    MPI_Bcast(&xPrime.front(), xPrime.size(), MPI_WEIGHT_T, 0, MPI_COMM_WORLD);
  }
#endif

  scrf->setParameters(xPrime);
  
  weight_t functionVal = 0;

  /* calculate our function value contribution */
  if (training) {
    for (int i=0; i<trainingAlignments.size(); i++) {
      myTrainingWork[i].startTime = time(NULL);
      
      if (objective == LIKELIHOOD)
	functionVal += computeLikelihood(i, true);
      else if (objective == EA)
	functionVal += computeEA(i, true);
      else if (objective == A)
	functionVal += computeA(i, true);
      else if (objective == PPP)
	functionVal += computePPP(i, true);
      else if (objective == ESS)
	functionVal += computeESS(i, true);
      else if (objective == ESA)
	functionVal += computeESA(i, true);
      else if (objective == ASA)
	functionVal += computeASA(i, true);
      
      myTrainingWork[i].endTime = time(NULL);
      myTrainingWork[i].elapsedTime += myTrainingWork[i].endTime - myTrainingWork[i].startTime;
    }
  }
  else {
    for (int i=0; i<testingAlignments.size(); i++) {
      myTestingWork[i].startTime = time(NULL);
      
      if (objective == LIKELIHOOD)
	functionVal += computeLikelihood(i, false);
      else if (objective == EA)
	functionVal += computeEA(i, false);
      else if (objective == A)
	functionVal += computeA(i, false);
      else if (objective == PPP)
	functionVal += computePPP(i, false);
      else if (objective == ESS)
	functionVal += computeESS(i, false);
      else if (objective == ESA)
	functionVal += computeESA(i, false);
      else if (objective == ASA)
	functionVal += computeASA(i, false);
      
      
      myTestingWork[i].endTime = time(NULL);
      myTestingWork[i].elapsedTime += myTestingWork[i].endTime - myTestingWork[i].startTime;
    }
  }
  
#ifdef MULTI
  if (id == 0) {
    /* receive function value contributions from other processes and add them to ours */
    weight_t combinedFunctionVals;
    MPI_Reduce(&functionVal, &combinedFunctionVals, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
    functionVal = combinedFunctionVals;

    if (training && stochastic) {
      /* receive length contributions from other processes */
      int totalLength = 0;
      for (int i = 0; i < trainingAlignments.size(); i++)
	totalLength += trainingAlignments[i]->length;
      
      MPI_Reduce(&totalLength, &trainingLength, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    }
  }
  else {    
    /* send function value contribution to the root process */
    MPI_Reduce(&functionVal, NULL, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (training && stochastic) {
      /* send the length of our current training example to the root process */
      int totalLength = 0;
      for (int i = 0; i < trainingAlignments.size(); i++)
	totalLength += trainingAlignments[i]->length;
      
      MPI_Reduce(&totalLength, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  }
#endif

  
  if (id == 0) {
    //normalize by length
    if (training) {
      functionVal /= trainingLength;
    }
    else {
      functionVal /= testingLength;
    }

    // add regularization term
    functionVal += regularizationTerm(x);

    cerr << "computeObjective() completed in " << time(NULL) - startTime << " seconds" << endl;
  }

  return functionVal;
}

weight_t Optimizer::computeLikelihood(int seqNumber, bool training) {
  if (training)
    scrf->setSequences(trainingAlignments[seqNumber], trainingNumericSeqs[seqNumber], trainingESTSeqs[seqNumber]);
  else
    scrf->setSequences(testingAlignments[seqNumber], testingNumericSeqs[seqNumber], testingESTSeqs[seqNumber]);

  if (training)
    scrf->precomputeForDP(&(trainingSegmentations[seqNumber]));
  else
    scrf->precomputeForDP(&(testingSegmentations[seqNumber]));
  weight_t trainingPartition = scrf->forward();
  
  scrf->dpMatrix.freeEverything();

  scrf->precomputeForDP(NULL);
  weight_t allPartition = scrf->forward();

  weight_t likelihood = trainingPartition - allPartition; /* minimize negative log-likelihood */

  return (weight_t)-1 * likelihood; /* minimize negative log likelihood */
}

weight_t Optimizer::computeEA(int seqNumber, bool training) {
  if (training)
    scrf->setSequences(trainingAlignments[seqNumber], trainingNumericSeqs[seqNumber], trainingESTSeqs[seqNumber]);
  else
    scrf->setSequences(testingAlignments[seqNumber], testingNumericSeqs[seqNumber], testingESTSeqs[seqNumber]);

  scrf->precomputeForDP(NULL);

#pragma omp parallel sections
 { 
   #pragma omp section
   scrf->forward();
  
   #pragma omp section
   scrf->backward();
 }

 vector<weight_t> posteriors;
 weight_t ea;
 if (training) {
   ea = scrf->computeAnnotationPosteriors(trainingSegmentations[seqNumber], posteriors, lambda);
 }
 else {
   ea = scrf->computeAnnotationPosteriors(testingSegmentations[seqNumber], posteriors, lambda);
 }

 return (weight_t)-1 * ea; /* minimize negative expected accuracy */
}

weight_t Optimizer::computeA(int seqNumber, bool training) {
  if (training)
    scrf->setSequences(trainingAlignments[seqNumber], trainingNumericSeqs[seqNumber], trainingESTSeqs[seqNumber]);
  else
    scrf->setSequences(testingAlignments[seqNumber], testingNumericSeqs[seqNumber], testingESTSeqs[seqNumber]);

  scrf->precomputeForDP(NULL);
  
#pragma omp parallel sections
 {
  #pragma omp section
    scrf->forward();
  #pragma omp section
    scrf->backward();
 }

  scrf->computePositionalPosteriors();

  Segmentation* segmentation;
  if (training)
    segmentation = &(trainingSegmentations[seqNumber]);
  else
    segmentation = &(testingSegmentations[seqNumber]);

  vector<weight_t> posteriorDifference;
  scrf->computePosteriorDifference(*segmentation, posteriorDifference);

  weight_t aa = 0;
  for (pos_t j=0; j<scrf->alignment->length; j++) {
    if (scrf->states[segmentation->labels[j]].optimize)
      aa += lambda * Q(posteriorDifference[j], gamma);
    else
      aa += Q(posteriorDifference[j], gamma);
  }

  if (aa > (weight_t)scrf->alignment->length) return 0;  //numerical error, return a bad value

  return (weight_t)-1 * aa; /* minimize negative approximate accuracy */
}

weight_t Optimizer::computePPP(int seqNumber, bool training) {
  if (training)
    scrf->setSequences(trainingAlignments[seqNumber], trainingNumericSeqs[seqNumber], trainingESTSeqs[seqNumber]);
  else
    scrf->setSequences(testingAlignments[seqNumber], testingNumericSeqs[seqNumber], testingESTSeqs[seqNumber]);

  scrf->precomputeForDP(NULL);

#pragma omp parallel sections
 { 
   #pragma omp section
   scrf->forward();
  
   #pragma omp section
   scrf->backward();
 }

 vector<weight_t> posteriors;
 if (training)
   scrf->computeAnnotationPosteriors(trainingSegmentations[seqNumber], posteriors, lambda);
 else
   scrf->computeAnnotationPosteriors(testingSegmentations[seqNumber], posteriors, lambda);
 weight_t ppp = 0;
 for (int i=0; i<posteriors.size(); i++)
   ppp += posteriors[i];

 return (weight_t)-1 * ppp; /* minimize negative sum of log positional posteriors */
}

weight_t Optimizer::computeESS(int seqNumber, bool training) {
  if (training)
    scrf->setSequences(trainingAlignments[seqNumber], trainingNumericSeqs[seqNumber], trainingESTSeqs[seqNumber]);
  else
    scrf->setSequences(testingAlignments[seqNumber], testingNumericSeqs[seqNumber], testingESTSeqs[seqNumber]);

  scrf->precomputeForDP(NULL);

#pragma omp parallel sections
 { 
   #pragma omp section
   scrf->forward();
  
   #pragma omp section
   scrf->backward();
 }

 weight_t ess;
 if (training) {
   ess = scrf->computeESS(trainingSegmentations[seqNumber], lambda);
 }
 else {
   ess = scrf->computeESS(testingSegmentations[seqNumber], lambda);
 }

 return (weight_t)-1 * ess; /* minimize negative expected accuracy */
}

weight_t Optimizer::computeESA(int seqNumber, bool training) {
  if (training)
    scrf->setSequences(trainingAlignments[seqNumber], trainingNumericSeqs[seqNumber], trainingESTSeqs[seqNumber]);
  else
    scrf->setSequences(testingAlignments[seqNumber], testingNumericSeqs[seqNumber], testingESTSeqs[seqNumber]);

  scrf->precomputeForDP(NULL);

#pragma omp parallel sections
 { 
   #pragma omp section
   scrf->forward();
  
   #pragma omp section
   scrf->backward();
 }

 weight_t esa;
 if (training) {
   esa = scrf->computeESA(trainingSegmentations[seqNumber], lambda);
 }
 else {
   esa = scrf->computeESA(testingSegmentations[seqNumber], lambda);
 }

 return esa;
}

weight_t Optimizer::computeASA(int seqNumber, bool training) {
  if (training)
    scrf->setSequences(trainingAlignments[seqNumber], trainingNumericSeqs[seqNumber], trainingESTSeqs[seqNumber]);
  else
    scrf->setSequences(testingAlignments[seqNumber], testingNumericSeqs[seqNumber], testingESTSeqs[seqNumber]);

  scrf->precomputeForDP(NULL);

#pragma omp parallel sections
  { 
#pragma omp section
    scrf->forward();
    
#pragma omp section
    scrf->backward();
  }

  weight_t asa;
  weight_t feMultiplier;

  if (training) {
    asa = scrf->computeASA(trainingSegmentations[seqNumber], lambda, gamma, feMultiplier);
  }
  else {
    asa = scrf->computeASA(testingSegmentations[seqNumber], lambda, gamma, feMultiplier);
  }
  
  return asa;
}

weight_t Optimizer::regularizationTerm(const vector<weight_t> &x) {
  weight_t reg = 0;

  for (int i=0; i<scrf->states.size(); i++)
    reg += C * scrf->states[i].startWeight * scrf->states[i].startWeight;

  for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
    SequenceFeatureSet* sfs = scrf->sequenceFeatureSets[i];
    if (sfs->type == DNAKMER || sfs->type == DNAKMERPAIR) {
      for (int k=0; k<sfs->paramToWeight.size(); k++) {
	reg += C * sfs->paramToWeight[k] * sfs->paramToWeight[k];
      }
    }
  }

  //cerr << "Regularization term is " << reg << ", regularization term divided by C is " << reg / C << endl;

  return reg;
}

void Optimizer::regularizationGradient(const vector<weight_t> &x, vector<weight_t> &g) {
  g.clear();
  g.resize(x.size(), 0);

  for (int i=0; i<scrf->states.size(); i++)
    g[scrf->states[i].startParamID] += 2.0 * C * scrf->states[i].startWeight;

  for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
    SequenceFeatureSet* sfs = scrf->sequenceFeatureSets[i];
    if (sfs->type == DNAKMER || sfs->type == DNAKMERPAIR) {
      for (int k=0; k<sfs->paramToWeight.size(); k++) {
	g[sfs->globalParamID[k]] += 2.0 * C * sfs->paramToWeight[k];
      }
    }
  }
}

void Optimizer::setWorkAssignment(vector<int>& assignment, bool training) {
  if (training) {
    for (int i=0; i<trainingAlignments.size(); i++)
      delete trainingAlignments[i];
    trainingAlignments.clear();

    for (int i=0; i<trainingNumericSeqs.size(); i++)
      delete trainingNumericSeqs[i];
    trainingNumericSeqs.clear();

    for (int i=0; i<trainingESTSeqs.size(); i++)
      delete trainingESTSeqs[i];
    trainingESTSeqs.clear();

    trainingGTFs.clear();
    trainingSegmentations.clear();
    myTrainingWork.clear();
  }
  else {
    for (int i=0; i<testingAlignments.size(); i++)
      delete testingAlignments[i];
    testingAlignments.clear();

    for (int i=0; i<testingNumericSeqs.size(); i++)
      delete testingNumericSeqs[i];
    testingNumericSeqs.clear();

    for (int i=0; i<testingESTSeqs.size(); i++)
      delete testingESTSeqs[i];
    testingESTSeqs.clear();

    testingGTFs.clear();
    testingSegmentations.clear();
    myTestingWork.clear();
  }

  //old sequences have been deleted
  scrf->alignment = NULL;
  scrf->numericSeq = NULL;
  scrf->estSeq = NULL;

  //determine which sequences we'll need to load
  set<string> seqNames;  
  scrf->getNeededSequences(seqNames);

  for (int i=0; i<assignment.size() && assignment[i] != -1; i++) {
    if (id == 0) {
      cerr << "Process " << id << " processing assigned ";
      if (training)
	cerr << "training";
      else
	cerr << "testing";
      cerr << " sequence " << (i+1) << " of " << assignment.size() << endl;
    }

    int assignedID = assignment[i];

    if (training) {
      trainingAlignments.push_back(new AlignmentSequence(trainingAlignmentFiles[assignedID], seqNames));
      trainingNumericSeqs.push_back(scrf->createNumericSequence(trainingAlignments.back()));
      if (trainingESTSeqFiles[assignedID] == "None")
	trainingESTSeqs.push_back(NULL);
      else
	trainingESTSeqs.push_back(new ESTSequence(trainingESTSeqFiles[assignedID]));

      scrf->setSequences(trainingAlignments.back(), trainingNumericSeqs.back(), trainingESTSeqs.back());
      GTF singleGTF(trainingGTFFiles[assignedID]);
      singleGTF.discardOverlaps();
      singleGTF.discardShortExons(3);  //minmum exon length for training data is 3
      GTF multipleGTF(trainingGTFFiles[assignedID]);

      //if our CRF does not support 5' or 3' UTRs, delete those annotations from the GTFs
      if (scrf->getStateID("FivePrimeUTR") == -1 && scrf->getStateID("FivePrimeUTRSingle") == -1) {
	singleGTF.removeFeatureType("5UTR");
	multipleGTF.removeFeatureType("5UTR");
      }
      if (scrf->getStateID("ThreePrimeUTR") == -1) {
	singleGTF.removeFeatureType("3UTR");
	multipleGTF.removeFeatureType("3UTR");
      }

      trainingGTFs.push_back(multipleGTF);
      trainingSegmentations.push_back(Segmentation(singleGTF, *scrf));
      if (! trainingSegmentations.back().setLabels())
	cerr << "Error on " << trainingGTFFiles[assignedID] << endl;
      myTrainingWork.push_back(WorkUnit(assignedID, trainingAlignmentFiles[assignedID],
					trainingGTFFiles[assignedID], id));
    }
    else {
      testingAlignments.push_back(new AlignmentSequence(testingAlignmentFiles[assignedID], seqNames));
      testingNumericSeqs.push_back(scrf->createNumericSequence(testingAlignments.back()));
      if (testingESTSeqFiles[assignedID] == "None")
	testingESTSeqs.push_back(NULL);
      else
	testingESTSeqs.push_back(new ESTSequence(testingESTSeqFiles[assignedID]));
      scrf->setSequences(testingAlignments.back(), testingNumericSeqs.back(), testingESTSeqs.back());
      GTF singleGTF(testingGTFFiles[assignedID]);
      singleGTF.discardOverlaps(); 
      singleGTF.discardShortExons(3);  //minimum exon length is 3
      GTF multipleGTF(testingGTFFiles[assignedID]);

      //if our CRF does not support 5' or 3' UTRs, delete those annotations from the GTFs
      if (scrf->getStateID("FivePrimeUTR") == -1 && scrf->getStateID("FivePrimeUTRSingle") == -1) {
	singleGTF.removeFeatureType("5UTR");
	multipleGTF.removeFeatureType("5UTR");
      }
      if (scrf->getStateID("ThreePrimeUTR") == -1) {
	singleGTF.removeFeatureType("3UTR");
	multipleGTF.removeFeatureType("3UTR");
      }

      testingGTFs.push_back(multipleGTF);
      testingSegmentations.push_back(Segmentation(singleGTF, *scrf));
      if (! testingSegmentations.back().setLabels())
	cerr << "Error on " << testingGTFFiles[assignedID] << endl;
      myTestingWork.push_back(WorkUnit(assignedID, testingAlignmentFiles[assignedID],
				       testingGTFFiles[assignedID], id));
    }
  }
}

void Optimizer::resetWorkAssignment() {
  vector<int> trainingAssignment;
  vector<int> testingAssignment;

  for (int i=0; i<myTrainingWork.size(); i++)
    trainingAssignment.push_back(myTrainingWork[i].id);
  for (int i=0; i<myTestingWork.size(); i++)
    testingAssignment.push_back(myTestingWork[i].id);

  setWorkAssignment(trainingAssignment, true);
  setWorkAssignment(testingAssignment, false);
}

void Optimizer::computeScalingFactors() {

  //the scaling factors are just the number of occurences of each feature in the training data
  scalingFactors.resize(scrf->parameterMap.size(), 0);

  for (int i=0; i<trainingAlignments.size(); i++) {
    vector<weight_t> trainingExpectations(scrf->parameterMap.size());

    scrf->setSequences(trainingAlignments[i], trainingNumericSeqs[i], trainingESTSeqs[i]);

    scrf->precomputeForDP(&(trainingSegmentations[i]));
    scrf->backward();
    scrf->forward();
    scrf->featureExpectations(trainingExpectations);
    
    scrf->dpMatrix.freeEverything();

    for (int i=0; i<scalingFactors.size(); i++) {
      scalingFactors[i] += trainingExpectations[i];
    }
  }

  if (id == 0) {
    /* add smoothing term */
    for (int i=0; i<scalingFactors.size(); i++) {
      scalingFactors[i] += 1;
    }
  }

#ifdef MULTI
  if (id == 0) {
    /* receive contributions from other processes and add them to ours */
    vector<weight_t> combinedScalingFactors(scalingFactors.size(), 0);
    MPI_Reduce(&scalingFactors.front(), &combinedScalingFactors.front(), scalingFactors.size(), 
	       MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
    scalingFactors = combinedScalingFactors;

    weight_t totalCounts = 0;
    for (int i=0; i<scalingFactors.size(); i++)
      totalCounts += scalingFactors[i];

    /* pool the scaling factors by parameter type */
    weight_t startWeightScale = 0;
    weight_t transitionScale = 0;
    weight_t sfsScale = 0;
    weight_t transitionSFSScale = 0;

    for (int i=0; i<scrf->states.size(); i++)
      startWeightScale += scalingFactors[scrf->states[i].startParamID];
    for (int i=0; i<scrf->transitions.size(); i++)
      transitionScale += scalingFactors[scrf->transitions[i].globalParamID];
    for (int i=0; i<scrf->stateTypes.size(); i++) {
      for (int j=0; j<scrf->stateTypes[i].sequenceFeatureSets.size(); j++) {
	SequenceFeatureSet* sfs = scrf->stateTypes[i].sequenceFeatureSets[j];
	for (int k=0; k<sfs->globalParamID.size(); k++)
	  sfsScale += scalingFactors[sfs->globalParamID[k]];
      }
    }
    for (int i=0; i<scrf->transitionTypes.size(); i++) {
      for (int j=0; j<scrf->transitionTypes[i].sequenceFeatureSets.size(); j++) {
	SequenceFeatureSet* sfs = scrf->transitionTypes[i].sequenceFeatureSets[j];
	for (int k=0; k<sfs->globalParamID.size(); k++)
	  transitionSFSScale += scalingFactors[sfs->globalParamID[k]];	
      }
    }

    startWeightScale /= totalCounts;
    transitionScale /= totalCounts;
    sfsScale /= totalCounts;
    transitionSFSScale /= totalCounts;

    for (int i=0; i<scrf->states.size(); i++)
      scalingFactors[scrf->states[i].startParamID] = startWeightScale;
    for (int i=0; i<scrf->transitions.size(); i++)
      scalingFactors[scrf->transitions[i].globalParamID] = transitionScale;
    for (int i=0; i<scrf->stateTypes.size(); i++) {
      for (int j=0; j<scrf->stateTypes[i].sequenceFeatureSets.size(); j++) {
	SequenceFeatureSet* sfs = scrf->stateTypes[i].sequenceFeatureSets[j];
	for (int k=0; k<sfs->globalParamID.size(); k++)
	  scalingFactors[sfs->globalParamID[k]] = sfsScale;
      }
    }
    for (int i=0; i<scrf->transitionTypes.size(); i++) {
      for (int j=0; j<scrf->transitionTypes[i].sequenceFeatureSets.size(); j++) {
	SequenceFeatureSet* sfs = scrf->transitionTypes[i].sequenceFeatureSets[j];
      for (int k=0; k<sfs->globalParamID.size(); k++)
	scalingFactors[sfs->globalParamID[k]] = transitionSFSScale;	
      }
    }
  }
  else {
    /* send contributions to root process */
    MPI_Reduce(&scalingFactors.front(), NULL, scalingFactors.size(), MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
  }
#endif
}

void Optimizer::scaleGridSearch(vector<weight_t>& initialGuess) {
  Scaler scaler(this, initialGuess);
  vector<weight_t> zero(1, 0.0);
  vector<weight_t> one(1, 1.0);
  weight_t scalingFactor = scaler.GoldenSection(zero, one, 0, 3, 10);
  for (int i=0; i<scrf->transitionTypes.size(); i++) {
    for (int k=0; k<scrf->transitionTypes[i].sequenceFeatureSets.size(); k++) {
      SequenceFeatureSet& sfs = *(scrf->transitionTypes[i].sequenceFeatureSets[k]);
      if (sfs.type == SVM) {
	for (int j=0; j<sfs.paramToWeight.size(); j++) {
	  initialGuess[sfs.globalParamID[j]] *= scalingFactor;
	}
	if (scrf->transitionTypes[i].name == "DonorSpliceGC") {
	  //take into account the fact that roughly 1% of donor sites are GC
	  for (int j=0; j<sfs.paramToWeight.size(); j++) {
	    initialGuess[sfs.globalParamID[j]] += log(0.01);
	  }
	}
      }
    }
  }
  return;


    weight_t stateScaleLower = 1.0;
    weight_t stateScaleUpper = 1.0;
    weight_t svmStateStep = 0.1;
    weight_t transitionScaleLower = 1;
    weight_t transitionScaleUpper = 5;
    weight_t svmTransitionStep = 0.5;
    
    weight_t minObjective = (weight_t)-1 * LOG_ZERO;  //infinity
    vector<weight_t> bestParameters;

    for (weight_t stateScale=stateScaleLower; stateScale<=stateScaleUpper; 
	 stateScale += svmStateStep) {
      for (weight_t transitionScale=transitionScaleLower; transitionScale<=transitionScaleUpper; 
	   transitionScale += svmTransitionStep) {

	vector<weight_t> testParameters = initialGuess;

	for (int i=0; i<scrf->stateTypes.size(); i++) {
	  for (int k=0; k<scrf->stateTypes[i].sequenceFeatureSets.size(); k++) {
	    SequenceFeatureSet& sfs = *(scrf->stateTypes[i].sequenceFeatureSets[k]);
	    if (sfs.type == SVM || sfs.type == DNAKMERPAIR) {
	      for (int j=0; j<sfs.paramToWeight.size(); j++) {
		testParameters[sfs.globalParamID[j]] *= stateScale;
	      }
	    }
	  }
	}
	for (int i=0; i<scrf->transitionTypes.size(); i++) {
	  for (int k=0; k<scrf->transitionTypes[i].sequenceFeatureSets.size(); k++) {
	    SequenceFeatureSet& sfs = *(scrf->transitionTypes[i].sequenceFeatureSets[k]);
	    if (sfs.type == SVM) {
	      for (int j=0; j<sfs.paramToWeight.size(); j++) {
		testParameters[sfs.globalParamID[j]] *= transitionScale;
	      }
	    }
	  }
	}

	weight_t testObjective = computeObjective(testParameters,true);
	cerr << "Scaling grid search (" << stateScale << "," << transitionScale << ") = " << testObjective << endl;
	if (testObjective < minObjective) {
	  minObjective = testObjective;
	  bestParameters = testParameters;

	  /*  GTFEvaluation trainingEval;
	  GTFEvaluation testingEval;
	  weight_t trainingLabelAccuracy;
	  weight_t testingLabelAccuracy;
	  predict(testParameters, trainingEval, testingEval, trainingLabelAccuracy, testingLabelAccuracy);
	  trainingLabelAccuracy /= trainingLength;
	  testingLabelAccuracy /= testingLength;
	  trainingLog << "GRID SEARCH (" << stateScale << "," << transitionScale << ") = " << testObjective << endl;
	  trainingLog << "Label accuracy on the training data is " << trainingLabelAccuracy << endl;
	  trainingLog << "Label accuracy on the testing data is " << testingLabelAccuracy << endl << endl;
	  
	  trainingLog << "ACCURACY ON TRAINING DATA" << endl;
	  trainingLog << trainingEval;
	  trainingLog << endl;
	  
	  trainingLog << "ACCURACY ON TESTING DATA" << endl;
	  trainingLog << testingEval;
	  trainingLog << endl << endl << endl << endl;
	  */
	}
      }
    }
    
    //set initial guess to the best set of parameters we found in our search
    initialGuess = bestParameters;
}

void Optimizer::printWorkStats(vector<int>& workSums, ostream& os) {
  weight_t minWork = workSums[0];
  weight_t maxWork = workSums[0];
  weight_t meanWork = 0;
  for (int i=0; i<workSums.size(); i++) {
    if ((weight_t)workSums[i] < minWork)
      minWork = workSums[i];
    if ((weight_t)workSums[i] > maxWork)
      maxWork = workSums[i];
    meanWork += workSums[i];
  }
  meanWork /= workSums.size();
  weight_t wastedPercent = (weight_t)100.0 * (maxWork - meanWork) / maxWork;
  os << setiosflags(ios::fixed | ios::showpoint) << setprecision(1);
  os << "Maximum work " << maxWork << ", minimum work " << minWork
     << ", mean work " << meanWork << endl;
  os << ((weight_t)100 - wastedPercent) << " percent load balancing efficiency" << endl;
  os << setprecision(7);
}


void Optimizer::makeInitialGuess(weight_t pseudocount, vector<weight_t>& initialGuess, bool newInitialGuess) {

  scrf->setParameterMap();

  initialGuess.clear();
  initialGuess.resize(scrf->parameterMap.size());

  /* set the initial guess based on a MLE with the assumption that
     all the features are independent (HMM-style) */
  vector<weight_t> numerators(initialGuess.size(), 0);
  vector<weight_t> denominators(initialGuess.size(), 0);
  
  for (int i=0; i<trainingAlignments.size(); i++) {
    scrf->setSequences(trainingAlignments[i], trainingNumericSeqs[i], trainingESTSeqs[i]);
    for (int j=0; j<trainingSegmentations[i].segments.size(); j++) {
      Segment& segment = trainingSegmentations[i].segments[j];
      State& state = scrf->states[segment.state];
      StateType& stateType = scrf->stateTypes[state.typeID];
      
      /* add counts for start parameters */
      weight_t length = segment.end - segment.start + 1;
      for (int k=0; k<scrf->states.size(); k++) {
	if (k == segment.state)
	  numerators[scrf->states[k].startParamID] += length;
	denominators[scrf->states[k].startParamID] += length;
      }
      if (segment.start != 0) {
	/* add counts for transition parameters */
	State& prevState = scrf->states[trainingSegmentations[i].segments[j-1].state];
	for (int k=0; k<prevState.transitionsFrom.size(); k++) {
	  TransitionFeature& transition = scrf->transitions[prevState.transitionsFrom[k]];
	  if (transition.toState == segment.state && 
	      transition.allowed(trainingAlignments[i], segment.start)) {
	    numerators[transition.globalParamID] += (weight_t)1;
	    if (transition.typeID != -1) {
	      /* add counts for transition sequence features */
	      TransitionType& transitionType = scrf->transitionTypes[transition.typeID];
	      for (int a=0; a<transitionType.sequenceFeatureSets.size(); a++) {
		SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[a];

		if (sfs->type == SVM) {
		  for (int b=0; b<sfs->valueToGlobalParamID.size(); b++) {
		    numerators[sfs->valueToGlobalParamID[b]] += (weight_t)1;
		    denominators[sfs->valueToGlobalParamID[b]] += (weight_t)1;
		  }
		}
		else {
		  sfs->positionExpectationHelper(segment.start, transition.strand, 1, numerators);
		  sfs->positionInitialGuessHelper(segment.start, transition.strand, denominators);
		}
	      }
	    }
	  }
	  denominators[transition.globalParamID] += (weight_t)1;
	}
      }
      /* add counts for all self transitions */
      for (pos_t pos=segment.start; pos<segment.end; pos++) {
	for (int k=0; k<state.transitionsFrom.size(); k++) {
	  TransitionFeature& transition = scrf->transitions[state.transitionsFrom[k]];
	  if (transition.toState == segment.state)
	    numerators[transition.globalParamID] += (weight_t)1;
	  denominators[transition.globalParamID] += (weight_t)1;
	}
      }
      
      for (int k=0; k<stateType.sequenceFeatureSets.size(); k++) {
	SequenceFeatureSet* sfs = stateType.sequenceFeatureSets[k];
	for (pos_t pos=segment.start; pos<=segment.end; pos++) {
	  sfs->positionExpectationHelper(pos, state.strand, 1, numerators);
	  sfs->positionInitialGuessHelper(pos, state.strand, denominators);
	}
      }
    }
  }

#ifdef MULTI
  if (id != 0) {
    /* send counts to the root process */
    MPI_Reduce(&numerators.front(), NULL, numerators.size(), MPI_WEIGHT_T,
	       MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&denominators.front(), NULL, denominators.size(), MPI_WEIGHT_T,
	       MPI_SUM, 0, MPI_COMM_WORLD);
  }
  else {
    /* collect counts from all other processes */
    vector<weight_t> myNumerators = numerators;
    vector<weight_t> myDenominators = denominators;
    MPI_Reduce(&myNumerators.front(), &numerators.front(), numerators.size(), 
	       MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myDenominators.front(), &denominators.front(), denominators.size(), 
	       MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
  }
#endif

  if (id == 0) {
    /* add pseudocounts */
    cerr << "Adding pseudocount of " << pseudocount << endl;
    for (int i=0; i<scrf->stateTypes.size(); i++) {  //lengths
      LengthFS* lfs = scrf->stateTypes[i].lengthFeatures;
      if (lfs != NULL && lfs->parameterization != "Geometric") {
	for (int j=1; j<=lfs->maxLength; j++) {
	  numerators[lfs->valueToGlobalParamID[j]] += pseudocount / (weight_t)(lfs->binSize);
	  denominators[lfs->valueToGlobalParamID[j]] += pseudocount / (weight_t)(lfs->binSize) * (weight_t)(lfs->maxLength);
	}
      }
    }
    for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) { //sequence feature sets
      SequenceFeatureSet* sfs = scrf->sequenceFeatureSets[i];
      for (int j=0; j<sfs->valueToGlobalParamID.size(); j++) {
	numerators[sfs->valueToGlobalParamID[j]] += pseudocount;
	set<int> denominatorsToIncrement;
	sfs->getDenominators(j, denominatorsToIncrement);
	for (set<int>::iterator iter = denominatorsToIncrement.begin(); iter != denominatorsToIncrement.end(); ++iter)
	  denominators[*iter] += pseudocount;
      }
    }

    /* coding state length distributions should normalize to 3 because the length % 3 is known
       from the state path */
    for (int i=0; i<scrf->stateTypes.size(); i++) {
      string typeName = scrf->stateTypes[i].name;
      if (typeName == "Initial" || typeName == "Internal" ||typeName == "Terminal" ||typeName == "Single") {
	LengthFS* lfs = scrf->stateTypes[i].lengthFeatures;
	if (lfs != NULL && lfs->parameterization != "Geometric") {
	  for (int j=0; j<lfs->globalParamID.size(); j++) {
	    numerators[lfs->globalParamID[j]] *= 3.0;
	  }
	}
      }
    }

    for (int i=0; i<initialGuess.size(); i++) {
      if (numerators[i] == (weight_t)0) {
	if (denominators[i] == (weight_t)0)
	  initialGuess[i] = (weight_t)0;
	else
	  initialGuess[i] = -30;  //very small
      } 
      else
	initialGuess[i] = log(numerators[i] / denominators[i]);
    }
    for (int i=0; i<scrf->stateTypes.size(); i++) {
      if (scrf->stateTypes[i].lengthFeatures != NULL &&
	  scrf->stateTypes[i].lengthFeatures->parameterization == "Geometric") {  /* special case */
	int baseParamID = scrf->stateTypes[i].lengthFeatures->valueToGlobalParamID[0];
	if (numerators[baseParamID] == (weight_t)0 || denominators[baseParamID] == (weight_t)0) {
	  /* State not observed, so transitions to it should be disallowed
	     Set these values to 0 anyway to avoid NaN errors */
	  initialGuess[baseParamID] = (weight_t)0;
	  initialGuess[baseParamID + 1] = (weight_t)0;
	}
	else {
	  initialGuess[baseParamID] = log(numerators[baseParamID] / denominators[baseParamID]);
	  initialGuess[baseParamID + 1] = log((weight_t)1 - numerators[baseParamID] / denominators[baseParamID]);
	}
      }
    }

    scrf->setParameters(initialGuess);
  }

  if (id == 0) {
    /* set SVM initial guess */
    for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
      if (scrf->sequenceFeatureSets[i]->type == SVM) {
	SVMFS* svmfs = (SVMFS*)scrf->sequenceFeatureSets[i];
	svmfs->makeInitialGuess();

	//delete all SVM training examples to save memory
	vector<vector<int> >().swap(svmfs->trainingTrue);
	vector<vector<int> >().swap(svmfs->trainingDecoys);
	vector<vector<int> >().swap(svmfs->testingTrue);
	vector<vector<int> >().swap(svmfs->testingDecoys);
      }
    }
  }

  if (id == 0) {
    //set some parameters to small random values
    for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
      SequenceFeatureSet* sfs = scrf->sequenceFeatureSets[i];
      if (sfs->type == MASKING || sfs->type == DNAKMER || sfs->type == DNAKMERPAIR || 
	  sfs->type == ESTPOSITION || sfs->type == ESTTRANSITION) {
	for (int j=0; j<sfs->paramToWeight.size(); j++) {
	  sfs->paramToWeight[j] = 1e-3 * rand() / RAND_MAX - 0.5e-3;
	}
      }
    }
  }


  scrf->setParameterMap();
  
  if (id == 0) {
    //check normalization for all distributions
    for (int i=0; i<scrf->stateTypes.size(); i++) {
      if (scrf->stateTypes[i].lengthFeatures != NULL) {
	weight_t sum = 0;
	for (int j=1; j<=scrf->stateTypes[i].lengthFeatures->maxLength; j++)
	  sum += exp(scrf->stateTypes[i].lengthFeatures->valueToWeight[j]);
	cerr << scrf->stateTypes[i].name << " length features normalize to " << sum << endl;
      }
    }
  }

  if (id == 0) {
    initialGuess.clear();
    initialGuess.resize(scrf->parameterMap.size());
    for (int i=0; i<initialGuess.size(); i++)
      initialGuess[i] = *(scrf->parameterMap[i]);
    
    cerr << "Made new initial guess" << endl;
  }
}

void Optimizer::testSVMScores() {

  for (int i=0; i<trainingAlignments.size(); i++) {
    scrf->setSequences(trainingAlignments[i], trainingNumericSeqs[i], trainingESTSeqs[i]);
    for (int j=0; j<trainingSegmentations[i].segments.size(); j++) {
      Segment& segment = trainingSegmentations[i].segments[j];
      State& state = scrf->states[segment.state];
      StateType& stateType = scrf->stateTypes[state.typeID];
      if (segment.start != 0) {
	State& prevState = scrf->states[trainingSegmentations[i].segments[j-1].state];
	for (int k=0; k<prevState.transitionsFrom.size(); k++) {
	  TransitionFeature& transition = scrf->transitions[prevState.transitionsFrom[k]];
	  if (transition.toState == segment.state && transition.typeID != -1) {
	    TransitionType& transitionType = scrf->transitionTypes[transition.typeID];
	    for (int a=0; a<transitionType.sequenceFeatureSets.size(); a++) {
	      SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[a];
	      if (sfs->type == SVM) {
		SVMFS* svmfs = (SVMFS*)sfs;
		
		cerr << scrf->states[transition.fromState].name << " --> "
		     << scrf->states[transition.toState].name << ", bin = ";
		if (transition.strand == PLUS || transition.strand == NONE)
		  cerr << scrf->numericSeq->tags[2 * svmfs->id][segment.start] << endl;
		else 
		  cerr << scrf->numericSeq->tags[2 * svmfs->id + 1][segment.start] << endl;
	      }
	    }
	  }
	}
      }
    }
  }
}

void Optimizer::getLengths() {
#ifdef MULTI
  if (id == 0) {
    //send command to all processes to compute lengths and send them to us
    crf_command_t command = SEND_LENGTHS;
    MPI_Bcast(&command, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif

  weight_t myTrainingLength = 0;
  for (int i=0; i<trainingAlignments.size(); i++)
    myTrainingLength += trainingAlignments[i]->length;
  
  weight_t myTestingLength = 0;
  for (int i=0; i<testingAlignments.size(); i++)
    myTestingLength += testingAlignments[i]->length;
  
  weight_t myTrainingPositionsToOptimize = 0;
  for (int i=0; i<trainingSegmentations.size(); i++) {
    for (int j=0; j<trainingSegmentations[i].labels.size(); j++) {
      if (scrf->states[trainingSegmentations[i].labels[j]].optimize)
	myTrainingPositionsToOptimize += 1.0;
    } 
  }

  weight_t myTestingPositionsToOptimize = 0;
  for (int i=0; i<testingSegmentations.size(); i++) {
    for (int j=0; j<testingSegmentations[i].labels.size(); j++) {
      if (scrf->states[testingSegmentations[i].labels[j]].optimize)
	myTestingPositionsToOptimize += 1.0;
    } 
  }

#ifdef MULTI
  if (id == 0) {
    trainingLength = 0;
    testingLength = 0;
    trainingPositionsToOptimize = 0;
    testingPositionsToOptimize = 0;
  }
  // collect all lengths at root process
  MPI_Reduce(&myTrainingLength, &trainingLength, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&myTestingLength, &testingLength, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&myTrainingPositionsToOptimize, &trainingPositionsToOptimize, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&myTestingPositionsToOptimize, &testingPositionsToOptimize, 1, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  trainingLength = myTrainingLength;
  testingLength = myTestingLength;
  trainingPositionsToOptimize = myTrainingPositionsToOptimize;
  testingPositionsToOptimize = myTestingPositionsToOptimize;
#endif
}	

#ifdef MULTI
void Optimizer::sendCommand(crf_command_t command, int destinationID) {
  MPI_Send(&command, 1, MPI_INT, destinationID, 0, MPI_COMM_WORLD);
}

void Optimizer::broadcastCommand(crf_command_t command) {
  MPI_Bcast(&command, 1, MPI_INT, 0, MPI_COMM_WORLD);  
}
#endif

/* checks to make sure forward and backward match up */
void Optimizer::testForwardBackward() {
  for (int i=0; i<trainingAlignments.size(); i++) {
    scrf->setSequences(trainingAlignments[i], trainingNumericSeqs[i], trainingESTSeqs[i]);

    if (objective == LIKELIHOOD) {
      weight_t forwardPartition;
      weight_t backwardPartition;  
      weight_t alphaPrimePartition;

      scrf->precomputeForDP(&(trainingSegmentations[i]));
      backwardPartition = scrf->backward();
      forwardPartition = scrf->forward();
      //alphaPrimePartition = scrf->computeAlphaPrime();
      
      cerr << "Training sequence " << i << " (fixed path), forward partition " << forwardPartition 
         << ", backward partition " << backwardPartition << endl; 
      //cerr << "Alpha prime partition " << alphaPrimePartition;
      scrf->dpMatrix.freeEverything();
      
      scrf->precomputeForDP(NULL);
      backwardPartition = scrf->backward();
      forwardPartition = scrf->forward();
      //alphaPrimePartition = scrf->computeAlphaPrime();      

      cerr << "Training sequence " << i << ", forward partition " << forwardPartition 
	   << ", backward partition " << backwardPartition << endl;
      //cerr << "Alpha prime partition " << alphaPrimePartition;
    }
    else if (objective == EA) {
      weight_t partition;
      weight_t alphaBetaEA;
      weight_t alphaStarEA;
      weight_t betaStarEA;
      vector<weight_t> posteriors;
      
      scrf->precomputeForDP(NULL);
      partition = scrf->forward();
      scrf->backward();
      scrf->computeAnnotationPosteriors(trainingSegmentations[i], posteriors, lambda);
      alphaBetaEA = 0.0;
      for (int j=0; j<posteriors.size(); j++) {
	if (scrf->states[trainingSegmentations[i].labels[j]].optimize)
	  alphaBetaEA += lambda * exp(posteriors[j]);
	else
	  alphaBetaEA += exp(posteriors[j]);
      }
      alphaStarEA = scrf->computeAlphaStar(trainingSegmentations[i], lambda);
      betaStarEA = scrf->computeBetaStar(trainingSegmentations[i], lambda);
      cerr << setiosflags(ios::fixed | ios::showpoint) << setprecision(30);
      cerr << "Partition = " << partition << endl;
      cerr << "Alpha/beta EA = " << alphaBetaEA << endl;
      cerr << "Alpha* EA = " << exp(alphaStarEA - partition) << endl;
      cerr << "Beta* EA = " << exp(betaStarEA - partition) << endl;
    }
    else if (objective == ESS) {
      weight_t partition;
      weight_t alphaBetaESS;
      weight_t alphaStarESS;
      weight_t betaStarESS;
      
      scrf->precomputeForDP(NULL);
      partition = scrf->forward();
      scrf->backward();
      alphaBetaESS = scrf->computeESS(trainingSegmentations[i], lambda);
      alphaStarESS = scrf->computeESSAlphaStar(trainingSegmentations[i], lambda);
      betaStarESS = scrf->computeESSBetaStar(trainingSegmentations[i], lambda);
      cerr << setiosflags(ios::fixed | ios::showpoint) << setprecision(30);
      cerr << "Partition = " << partition << endl;
      cerr << "Alpha/beta ESS = " << alphaBetaESS << endl;
      cerr << "Alpha* ESS = " << exp(alphaStarESS - partition) << endl;
      cerr << "Beta* ESS = " << exp(betaStarESS - partition) << endl;
    }
    else if (objective == ESA) {
      weight_t partition;
      weight_t alphaBetaESA;
      weight_t alphaStarESA;
      weight_t betaStarESA;
      
      scrf->precomputeForDP(NULL);
      partition = scrf->forward();
      scrf->backward();
      alphaBetaESA = scrf->computeESA(trainingSegmentations[i], lambda);

      scrf->dpMatrix.freeEverything();
      scrf->precomputeForDP(NULL);
      alphaStarESA = scrf->computeESAAlphaStar(trainingSegmentations[i], lambda);
      betaStarESA = scrf->computeESABetaStar(trainingSegmentations[i], lambda);

      cerr << setiosflags(ios::fixed | ios::showpoint) << setprecision(30);
      cerr << "Partition = " << partition << endl;
      cerr << "Alpha/beta ESA = " << alphaBetaESA << endl;
      cerr << "Alpha* ESA = " << exp(alphaStarESA - partition) << endl;
      cerr << "Beta* ESA = " << exp(betaStarESA - partition) << endl;
    }
    else if (objective == ASA) {
      weight_t partition;
      weight_t alphaBetaASA;
      weight_t alphaStarASA;
      weight_t betaStarASA;
      weight_t feMultiplier;

      scrf->precomputeForDP(NULL);
      partition = scrf->forward();
      scrf->backward();
      alphaBetaASA = scrf->computeASA(trainingSegmentations[i], lambda, gamma, feMultiplier);
      alphaStarASA = scrf->computeASAAlphaStar(trainingSegmentations[i], lambda, gamma);
      betaStarASA = scrf->computeASABetaStar(trainingSegmentations[i], lambda, gamma);

      cerr << setiosflags(ios::fixed | ios::showpoint) << setprecision(30);
      cerr << "Partition = " << partition << endl;
      cerr << "Alpha/beta ASA = " << alphaBetaASA << endl;
      cerr << "Alpha* ASA = " << exp(alphaStarASA - partition) << endl;
      cerr << "Beta* ASA = " << exp(betaStarASA - partition) << endl;
    }
  }
}

void Optimizer::testGradient() {

  trainingLog << "TESTING GRADIENT..." << endl;
  
  /* get the current parameters for the sCRF */
  vector<weight_t> params;
  for (int i=0; i<scrf->parameterMap.size(); i++)
    params.push_back(*(scrf->parameterMap[i]));

  vector<weight_t> gradient;

  weight_t ll = ComputeFunction(params);
  trainingLog << "Objective: " << ll << endl;
  trainingLog << "Objective from ComputeFunctionAndGradient: " << ComputeFunctionAndGradient(gradient, params) << endl;

  ComputeGradient(gradient, params);

  /*
  trainingLog << "GRADIENT: " << endl;
  for (int i=0; i<gradient.size(); i++) {
    trainingLog << i << "\t" << gradient[i];
    if (scaling)
      trainingLog << "\t" << scalingFactors[i] << endl;
    trainingLog << endl;
  } 
  */
  
  for (int i=0; i<gradient.size(); i++) {
    trainingLog << "Parameter " << i << endl;
    if (scaling)
      trainingLog << "Scaling factor: " << scalingFactors[i] << endl;
    bool computeFD = true;

    /* figure out parameter name */
    for (int j=0; j<scrf->states.size(); j++) {
      if (scrf->states[j].startParamID == i) {
	trainingLog << "State " << scrf->states[i].name << " start weight" << endl;
        computeFD = false;
      }
    }
    for (int j=0; j<scrf->transitions.size(); j++) {
      if (scrf->transitions[j].globalParamID == i) {
	trainingLog << "Transition from " << scrf->states[scrf->transitions[j].fromState].name <<
	  "\t to " << scrf->states[scrf->transitions[j].toState].name << endl;
	computeFD = false;
      }	
    }
    for (int j=0; j<scrf->stateTypes.size(); j++) {
      if (scrf->stateTypes[j].lengthFeatures != NULL) {
	for (int k=0; k<scrf->stateTypes[j].lengthFeatures->globalParamID.size(); k++) {
	  if (scrf->stateTypes[j].lengthFeatures->globalParamID[k] == i) {
	    trainingLog << "State type " << scrf->stateTypes[j].name << " length feature " << k << endl;
	    //computeFD = false;
	  }
	}
      }
    }

    for (int j=0; j<scrf->stateTypes.size(); j++) {
      for (int k=0; k<scrf->stateTypes[j].sequenceFeatureSets.size(); k++) {
	SequenceFeatureSet *sfs = scrf->stateTypes[j].sequenceFeatureSets[k];

	for (int a=0; a<sfs->paramToWeight.size(); a++) {
	  if (sfs->globalParamID[a] == i) {
	    trainingLog << "State type " << scrf->stateTypes[j].name << " sequenceFeatureSet type " <<
	      sfs->type << " parameter " << a << endl;
	  }
	}
      }
    }
    for (int j=0; j<scrf->transitionTypes.size(); j++) {
      for (int k=0; k<scrf->transitionTypes[j].sequenceFeatureSets.size(); k++) {
	SequenceFeatureSet *sfs = scrf->transitionTypes[j].sequenceFeatureSets[k];
	for (int a=0; a<sfs->paramToWeight.size(); a++) {
	  if (sfs->globalParamID[a] == i)
	    trainingLog << "Transition type " << scrf->transitionTypes[j].name << " sequenceFeatureSet type " <<
	      sfs->type << " parameter " << a << endl;
	  //computeFD = false;
	}
      }
    }

    weight_t originalParam = params[i];
    trainingLog << "Parameter value: " << originalParam << endl;
    if (computeFD) {
      cerr << "Parameter " << i << endl;
      trainingLog << "Gradient value: " << gradient[i] << endl;
      for (int j=0; j<=7; j++) {
	weight_t eps = 1;
	for (int k=0; k<j; k++)
	  eps *= 0.1;
	
	params[i] = originalParam + eps;
	weight_t new_ll = ComputeFunction(params);
	weight_t pd = (new_ll - ll) / eps;
	if (pd == (weight_t)0) {
	  trainingLog << "**" << gradient[i] << "/" << pd << "\t";
	}
	else {
	  weight_t relativeError = (gradient[i] - pd) / pd;
	  trainingLog << eps << "\t" << new_ll << "\t" << pd << " ("<< relativeError << ")\t";
	}
	trainingLog << endl;
      }
      trainingLog << endl << endl;
    }
    params[i] = originalParam;
  }
}


