#include "SVMTrainer.h"
#include <pthread.h>
#include <signal.h>

#ifdef MULTI
#include <mpi.h>
#endif

SVMTrainer::~SVMTrainer() {
  for (int i=0; i<trainingAlignments.size(); i++)
    delete trainingAlignments[i];
  for (int i=0; i<testingAlignments.size(); i++)
    delete testingAlignments[i];
}

SVMTrainer::SVMTrainer(SemiCRF* scrf, string trainingList, string testingList) {
  this->scrf = scrf;

  id = 0;
  numProcs = 1;

#ifdef MULTI
  MPI_Comm_rank(MPI_COMM_WORLD, &id);  /* get our id */
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs); /* get number of process */
#endif

  cerr << "Process " << id << " now beginning SVM training" << endl;

  /* parse the training list */
  ifstream trainingListStream(trainingList.c_str());
  if (! trainingListStream.is_open())
    fatalError("Could not open training list");
  while (! trainingListStream.eof()) {
    string s;
    trainingListStream >> s;
    if (s.length() > 0)
      trainingAlignmentFiles.push_back(s);
    trainingListStream >> s;  //discard EST seq file
    trainingListStream >> s;
    if (s.length() > 0)
      trainingGTFFiles.push_back(s);
  } 
  if (trainingAlignmentFiles.size() < numProcs)
    fatalError("Optimizer: More processes requested than training examples found");

  /* Assign training examples */
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
      testingListStream >> s; //discard EST seq file
      testingListStream >> s;
      if (s.length() > 0)
	testingGTFFiles.push_back(s);
    }
  }

  /* Assign testing examples */
  vector<int> testingAssignment;
  for (int i=id; i<testingAlignmentFiles.size(); i += numProcs)
    testingAssignment.push_back(i);
  setWorkAssignment(testingAssignment, false);
}

void SVMTrainer::train() {
  //cerr << "Process " << id << " collecting training examples" << endl;

  //collect training examples
  for (int i=0; i<trainingAlignments.size(); i++)
    collectExamples(trainingAlignments[i], trainingSegmentations[i], true);

  //collect testing examples
  for (int i=0; i<testingAlignments.size(); i++)
    collectExamples(testingAlignments[i], testingSegmentations[i], false);

  //cerr << "Process " << id << " done collecting examples" << endl;

#ifdef MULTI
  if (id == 0) {
    //receive SVM feature vectors from all other processes
    for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
      SequenceFeatureSet* sfs = scrf->sequenceFeatureSets[i];
      if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	for (int otherID=1; otherID<numProcs; otherID++) {
	  vector<vector<int> > otherTrainingTrue;
	  vector<vector<int> > otherTrainingDecoys;
	  vector<vector<int> > otherTestingTrue;
	  vector<vector<int> > otherTestingDecoys;
	  receiveSparseVectors(otherID, otherTrainingTrue);
	  receiveSparseVectors(otherID, otherTrainingDecoys);
	  receiveSparseVectors(otherID, otherTestingTrue);
	  receiveSparseVectors(otherID, otherTestingDecoys);
	  svmfs->trainingTrue.insert(svmfs->trainingTrue.end(), otherTrainingTrue.begin(), otherTrainingTrue.end()); 
	  svmfs->trainingDecoys.insert(svmfs->trainingDecoys.end(), otherTrainingDecoys.begin(), otherTrainingDecoys.end()); 
	  svmfs->testingTrue.insert(svmfs->testingTrue.end(), otherTestingTrue.begin(), otherTestingTrue.end()); 
	  svmfs->testingDecoys.insert(svmfs->testingDecoys.end(), otherTestingDecoys.begin(), otherTestingDecoys.end()); 
	}
	cerr << svmfs->trainingTrue.size() << " true examples, " << svmfs->trainingDecoys.size() << " decoys" << endl; 
	
	/*
	cerr << "True examples:" << endl;
	for (int j=0; j<10; j++) {
	  svmfs->printAlignmentBlock(cerr, svmfs->trainingTrue[j]);
	  cerr << endl;
	}
	cerr << "Decoy examples:" << endl;
	for (int j=0; j<10; j++) {
	  svmfs->printAlignmentBlock(cerr, svmfs->trainingDecoys[j]);
	  cerr << endl;
	}
	*/
      }
    }

    //orchestrate SVM training
    const weight_t gridLeft = -15;
    const weight_t gridRight = 0;
    const weight_t gridStep = 0.1;
    const int numGridPoints = (int)((gridRight - gridLeft)/gridStep);
    vector<weight_t> gridPoints(numGridPoints);
    for (int i=0; i<gridPoints.size(); i++)
      gridPoints[i] = exp(gridLeft + i*gridStep);

    for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
      SequenceFeatureSet* sfs = scrf->sequenceFeatureSets[i];
      if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;

	time_t startTime = time(NULL);	
	cerr << "Training SVM ID " << svmfs->id << endl;
	cerr << svmfs->trainingTrue.size() << " true examples, " 
	     << svmfs->trainingDecoys.size() << " decoys" << endl;
	
	//send training examples to all processes
	for (int otherID=1; otherID<numProcs; otherID++) {
	  sendCommand(RECEIVE_SVM_TRAINING_DATA, otherID);
	  MPI_Send(&(svmfs->id), 1, MPI_INT, otherID, 0, MPI_COMM_WORLD);
	}
	sendSparseVectors(-1, svmfs->trainingTrue);
	sendSparseVectors(-1, svmfs->trainingDecoys);
	sendSparseVectors(-1, svmfs->testingTrue);
	sendSparseVectors(-1, svmfs->testingDecoys);
	
	MPI_Status status;
	weight_t cvAccuracy;
	int gridPoint;
	weight_t maxCVAccuracy = 0;
	int maxGridPoint;
	vector<weight_t> accuracies(numGridPoints, -1);
	
	//send initial work
	for (int otherID=1; otherID<numProcs && otherID <= numGridPoints; otherID++) {
	  sendCommand(CV_SVM, otherID);
	  int initialPoint = otherID - 1;
	  MPI_Send(&initialPoint, 1, MPI_INT, otherID, 0, MPI_COMM_WORLD);
	  MPI_Send(&(gridPoints[initialPoint]), 1, MPI_WEIGHT_T, otherID, 0, MPI_COMM_WORLD);
	}
	
	//receive results and send more work if necessary
	int totalWork = numGridPoints;
	int workReceived = 0;
	int nextGridPoint = numProcs - 1;
	svm_command_t response;
	while (workReceived < totalWork) {
	  MPI_Recv(&response, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

	  if (response == CV_SVM) {
	    MPI_Recv(&gridPoint, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
	    MPI_Recv(&cvAccuracy, 1, MPI_WEIGHT_T, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);

	    cerr << "Process " << status.MPI_SOURCE << " sent results: " 
		 << "grid point " << gridPoint << ", C = " << gridPoints[gridPoint] 
		 << ", accuracy " << cvAccuracy << endl;

	    workReceived++;	  

	    //record result
	    accuracies[gridPoint] = cvAccuracy;
	  
	    if (cvAccuracy > maxCVAccuracy) {
	      //these parameters are the best so far
	      cerr << "***** New maximum CV accuracy: " << cvAccuracy << " at grid point " << gridPoint << endl; 
	      maxCVAccuracy = cvAccuracy;
	      maxGridPoint = gridPoint;
	      //get SVM parameters from sender
	      sendCommand(SEND_SVM, status.MPI_SOURCE);
	      svmfs->receiveSVM(status.MPI_SOURCE);
	    }
	    if (nextGridPoint < numGridPoints) {
	      //tell sender to test next grid point
	      sendCommand(CV_SVM, status.MPI_SOURCE);
	      MPI_Send(&nextGridPoint, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	      weight_t svmC = gridPoints[nextGridPoint];
	      MPI_Send(&svmC, 1, MPI_WEIGHT_T, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	      nextGridPoint++;
	    }

	    //if we've processed three points to either side of the maximum, declare victory
	    //require that those points have significantly different accuracies
	    weight_t tol = 1e-6;
	    bool maxFound = true;
	    for (int i=maxGridPoint+1; i<gridPoints.size() && i<=maxGridPoint+3; i++) {
	      if (accuracies[i] < 0 || fabs(accuracies[i] - maxCVAccuracy) < tol)
		maxFound = false;
	    }
	    for (int i=maxGridPoint-1; i>=0 && i>=maxGridPoint-3; i--) {
	      if (accuracies[i] < 0 || fabs(accuracies[i] - maxCVAccuracy) < tol)
		maxFound = false;
	    }
	    if (maxFound) {
	      //abort SVM training
	      cerr << "Maximum found, aborting SVM training" << endl;
	      svm_command_t abortCommand = ABORT_SVM_TRAINING;
	      for (int otherID=1; otherID<numProcs; otherID++)
		MPI_Send(&abortCommand, 1, MPI_INT, otherID, 0, MPI_COMM_WORLD);

	      //print out summary of results
	      cerr << "Grid search results:" << endl;
	      for (int i=0; i<gridPoints.size(); i++) {
		if (accuracies[i] > 0) {
		  cerr << i << "\t" << gridPoints[i] << "\t" << accuracies[i];
		  if (i == maxGridPoint)
		    cerr << "\t*****";
		  cerr << endl;
		}
	      }

	      //wait for acknowledgement from all other processes
	      vector<bool> ackReceived(numProcs, false);
	      int acks = 0;
	      while (acks < numProcs-1) {
		MPI_Recv(&abortCommand, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		if (abortCommand == ABORT_SVM_TRAINING) {
		  if (!ackReceived[status.MPI_SOURCE]) {
		    acks++;
		    ackReceived[status.MPI_SOURCE] = true;
		  }
		}
		else if (abortCommand == CV_SVM) {
		  MPI_Recv(&gridPoint, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
		  MPI_Recv(&cvAccuracy, 1, MPI_WEIGHT_T, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
		  abortCommand = ABORT_SVM_TRAINING;
		  MPI_Send(&abortCommand, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}
		else {
		  fatalError("Unknown command during abort ack phase of SVM training");
		}
	      }

	      break;
	    }
	  }
	  else if (response == ABORT_SVM_TRAINING) {
	    //extra ack, do nothing
	  }
	  else {
	    cerr << "Unexpected message: " << response << endl;
	    fatalError("Invalid response during SVM training");
	  }
	}
	
	//send the best SVM parameters to all other processes
	for (int otherID=1; otherID<numProcs; otherID++) {
	  sendCommand(RECEIVE_SVM, otherID);
	  MPI_Send(&(svmfs->id), 1, MPI_INT, otherID, 0, MPI_COMM_WORLD);             
	  svmfs->sendSVM(otherID);
	}
	cerr << "SVM training took " << (time(NULL) - startTime) << " seconds" << endl;
      }
    }
    
    cerr << "SVM training complete" << endl;
    
    //inform all other processes SVM training is complete
    for (int otherID=1; otherID<numProcs; otherID++)
      sendCommand(SVM_TRAINING_COMPLETE, otherID);
  }
  else {
    //send all SVM feature vectors to the root process
    for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
      SequenceFeatureSet* sfs = scrf->sequenceFeatureSets[i];
      if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	sendSparseVectors(0, svmfs->trainingTrue);
	sendSparseVectors(0, svmfs->trainingDecoys);
	sendSparseVectors(0, svmfs->testingTrue);
	sendSparseVectors(0, svmfs->testingDecoys);
      }
    }

    //wait for instructions to do SVM training
    SVMFS* svmfs;
    int gridPoint;
    weight_t svmC;
    while (true) {
      weight_t cvAccuracy;
      svm_command_t command;      
      MPI_Recv(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);  /* receive command */
      
      if (command == RECEIVE_SVM_TRAINING_DATA) {
	//receive svmID
	int svmID;
	MPI_Recv(&svmID, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
	svmfs = scrf->getSVMFSByID(svmID);

	receiveSparseVectors(-1, svmfs->trainingTrue);
	receiveSparseVectors(-1, svmfs->trainingDecoys);
	receiveSparseVectors(-1, svmfs->testingTrue);
	receiveSparseVectors(-1, svmfs->testingDecoys);
      }
      else if (command == CV_SVM) {
	//receive which grid point we are processing
	MPI_Recv(&gridPoint, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
	//receive parameter to use for training
	MPI_Recv(&svmC, 1, MPI_WEIGHT_T, 0, 0, MPI_COMM_WORLD, NULL);
	svmfs->C = svmC;

	//cerr << "Process " << id << " handling grid point " << gridPoint << ", training with C = " << svmC << endl;

	svmfs->cvSVM();

	//send results back to master process
	MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Send(&gridPoint, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Send(&(svmfs->cvAccuracy), 1, MPI_WEIGHT_T, 0, 0, MPI_COMM_WORLD);
      }
      else if (command == SEND_SVM) {
	svmfs->sendSVM(0);
      }
      else if (command == RECEIVE_SVM) {
	//receive svmID
	int svmID;
	MPI_Recv(&svmID, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
	svmfs = scrf->getSVMFSByID(svmID);
	svmfs->receiveSVM(0);
      }
      else if (command == SVM_TRAINING_COMPLETE) {
	//move on to CRF training
	cerr << "Process " << id << " finished SVM training" << endl;
	break;
      }
      else if (command == ABORT_SVM_TRAINING) {
	MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);  //send back abort command as ack
      }
      else {
	fatalError("Received invalid command during SVM training");
      }
    }
  }
#else
  //single processor mode, just train the SVMs with C = 1
  for (int i=0; i<scrf->sequenceFeatureSets.size(); i++) {
    SequenceFeatureSet* sfs = scrf->sequenceFeatureSets[i];
    if (sfs->type == SVM) {
      SVMFS* svmfs = (SVMFS*)sfs;
      svmfs->C = 1;
      svmfs->cvSVM();
    }
  }
#endif
}

//collects an example from noncoding sequence to use as a decoy for the given SVM
void SVMTrainer::addNoncodingDecoy(SVMFS* svmfs, AlignmentSequence* alignment, Segmentation& segmentation,
				   bool training) {
  while (true) { 
    //pick a position at random
    pos_t samplePos = (pos_t)((double)rand() / RAND_MAX * alignment->length);

    //if it falls in noncoding sequence, add a decoy and return
    for (pos_t j=0; j<segmentation.segments.size(); j++) {
      Segment& segment = segmentation.segments[j];
      if (samplePos >= segment.start && samplePos <= segment.end) {
	StateType& stateType = scrf->stateTypes[scrf->states[segment.state].typeID];
	if (stateType.name == "Noncoding") {
	  vector<int> x;
	  if (svmfs->buildInputVector(x, samplePos, scrf->states[segment.state].strand)) {
	    if (training)
	      svmfs->trainingDecoys.push_back(x);
	    else
	      svmfs->testingDecoys.push_back(x);
	    return;
	  }
	}
      }
    }
  }
}

//search for a decoy from the alignment matching the k-mer
//requirements of the given transition, starting close to the transition's position
void SVMTrainer::addCloseDecoy(SVMFS* svmfs, AlignmentSequence* alignment, 
			  int transitionID, pos_t pos, bool training) {
  //randomly select a search direction
  pos_t start;
  pos_t increment;
  if (2 * rand() / RAND_MAX == 0) {
    start = pos - 1;
    increment = -1;
  }
  else {
    start = pos + 1;
    increment = 1;
  }

  vector<int> x;
  bool success = false;
  pos_t j;
  for (j=start; j>=scrf->allowedKmerOffset && 
	 j<alignment->length - scrf->allowedKmerLength + scrf->allowedKmerOffset - 1; j += increment) {
    int kmer = 0;
    for (int kmerPos=0; kmerPos<scrf->allowedKmerLength; kmerPos++)
      kmer += DNA_INDEX_COEFF[scrf->allowedKmerLength - kmerPos - 1] *
	DNA_TO_INDEX[alignment->sequenceArray[0][j - scrf->allowedKmerOffset + kmerPos]];
    
    if (scrf->allowedLookup[kmer][transitionID]) {
      if (svmfs->buildInputVector(x, j, scrf->transitions[transitionID].strand)) {
	if (training)
	  svmfs->trainingDecoys.push_back(x);
	else
	  svmfs->testingDecoys.push_back(x);
	success = true;
	break;
      }
      //else
      //cerr << pos << ": transition allowed at position " << j << ", but buildInputVector returned false" << endl;
    }
  }
  //if (! success)
  //cerr << pos << ": failed to find decoy, final position was " << j << endl;
}

//collects true and decoy examples for each SVM from an alignment and segmentation
void SVMTrainer::collectExamples(AlignmentSequence* alignment, Segmentation& segmentation,
				 bool training) {

  scrf->setSequences(alignment, NULL, NULL);
  for (int j=0; j<segmentation.segments.size(); j++) {
    Segment& segment = segmentation.segments[j];
    State& state = scrf->states[segment.state];
    StateType& stateType = scrf->stateTypes[state.typeID];

    //state sequence features
    for (int a=0; a<stateType.sequenceFeatureSets.size(); a++) {
      SequenceFeatureSet* sfs = stateType.sequenceFeatureSets[a];
      if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	vector<int> x;
	for (pos_t pos=segment.start; pos<=segment.end; pos++) {
	  if ((weight_t)rand() / (weight_t)RAND_MAX < svmfs->sampleRate &&
	      svmfs->buildInputVector(x, pos, state.strand)) {
	    if (training)
	      svmfs->trainingTrue.push_back(x);
	    else
	      svmfs->testingTrue.push_back(x);
	    addNoncodingDecoy(svmfs, alignment, segmentation, training);
	  }
	}
      }
    }

    //transition sequence features
    if (segment.start != 0) {
      State& prevState = scrf->states[segmentation.segments[j-1].state];
      for (int k=0; k<prevState.transitionsFrom.size(); k++) {
	TransitionFeature& transition = scrf->transitions[prevState.transitionsFrom[k]];
	if (transition.toState == segment.state && transition.typeID != -1 &&
	    transition.allowed(alignment, segment.start)) {
	  TransitionType& transitionType = scrf->transitionTypes[transition.typeID];
	  for (int a=0; a<transitionType.sequenceFeatureSets.size(); a++) {
	    SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[a];
	    if (sfs->type == SVM) {
	      SVMFS* svmfs = (SVMFS*)sfs;
	      vector<int> x;
	      if ((weight_t)rand() / (weight_t)RAND_MAX < svmfs->sampleRate &&
		    svmfs->buildInputVector(x, segment.start, transition.strand)) {
		if (training)
		  svmfs->trainingTrue.push_back(x);
		else
		  svmfs->testingTrue.push_back(x);
		addCloseDecoy(svmfs, alignment, prevState.transitionsFrom[k], segment.start, training);
	      }
	    }
	  }
	}
      }
    }
  }
}


#ifdef MULTI
void SVMTrainer::sendCommand(svm_command_t command, int destinationID) {
  MPI_Send(&command, 1, MPI_INT, destinationID, 0, MPI_COMM_WORLD);
}

void SVMTrainer::broadcastCommand(svm_command_t command) {
  MPI_Bcast(&command, 1, MPI_INT, 0, MPI_COMM_WORLD);  
}

void SVMTrainer::sendSparseVectors(int destinationID, vector<vector<int> >& sparseVectors) {

  //calculate lenght of messages and inform destination(s)
  int length = 0;
  for (int i=0; i<sparseVectors.size(); i++)
    length += sparseVectors[i].size() + 1;
  if (destinationID == -1)  //broadcast
    MPI_Bcast(&length, 1, MPI_INT, id, MPI_COMM_WORLD);
  else
    MPI_Send(&length, 1, MPI_INT, destinationID, 0, MPI_COMM_WORLD);

  int* indices = (int*) malloc (length * sizeof(int));
  weight_t* values = (weight_t*) malloc (length * sizeof(weight_t));
  int k = 0;
  for (int i=0; i<sparseVectors.size(); i++) {
    for (int j=0; j<sparseVectors[i].size(); j++) {
      indices[k] = sparseVectors[i][j];
      values[k] = 1;
      k++;
    }
    indices[k] = -1;  //marks vector boundary
    values[k] = -1;
    k++;
  }
  if (destinationID == -1) {
    //broadcast
    MPI_Bcast(indices, length, MPI_INT, id, MPI_COMM_WORLD);
    MPI_Bcast(values, length, MPI_WEIGHT_T, id, MPI_COMM_WORLD);
  }
  else {
    MPI_Send(indices, length, MPI_INT, destinationID, 0, MPI_COMM_WORLD);
    MPI_Send(values, length, MPI_WEIGHT_T, destinationID, 0, MPI_COMM_WORLD);
  }

  free(indices);
  free(values);

  //cerr << "Process " << id << " sent " << sparseVectors.size() << " sparse vectors" << endl;
}

void SVMTrainer::receiveSparseVectors(int sourceID, vector<vector<int> >& sparseVectors) {
  
  int length;
  if (sourceID == -1)  //broadcast from the root process
    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
  else
    MPI_Recv(&length, 1, MPI_INT, sourceID, 0, MPI_COMM_WORLD, NULL);

  int* indices = (int*) malloc (length * sizeof(int));
  weight_t* values = (weight_t*) malloc (length * sizeof(weight_t));
  if (sourceID == -1) {
    //broadcast from root process
    MPI_Bcast(indices, length, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(values, length, MPI_WEIGHT_T, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Recv(indices, length, MPI_INT, sourceID, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(values, length, MPI_WEIGHT_T, sourceID, 0, MPI_COMM_WORLD, NULL);
  }

  int numVectors = 0;
  for (int i=0; i<length; i++) {
    if (indices[i] == -1)
      numVectors++;
  }

  sparseVectors.clear();
  sparseVectors.resize(numVectors);
  int k = 0;
  for (int i=0; i<length; i++) {
    if (indices[i] == -1)
      k++;
    else
      sparseVectors[k].push_back(indices[i]);
  }

  free(indices);
  free(values);

  //cerr << "Process " << id << " received " << sparseVectors.size() << " sparse vectors" << endl;
}
#endif


void SVMTrainer::setWorkAssignment(vector<int> assignment, bool training) {
  if (training) {
    for (int i=0; i<trainingAlignments.size(); i++)
      delete trainingAlignments[i];
    trainingAlignments.clear();
    trainingSegmentations.clear();
  }
  else {
    for (int i=0; i<testingAlignments.size(); i++)
      delete testingAlignments[i];
    testingAlignments.clear();
    testingSegmentations.clear();
  }

  scrf->alignment = NULL;  //old alignment has been deleted

  //figure out what sequences we'll need to load
  set<string> seqNames;  
  scrf->getNeededSequences(seqNames);

  for (int i=0; i<assignment.size() && assignment[i] != -1; i++) {
    int assignedID = assignment[i];

    if (training) {
      trainingAlignments.push_back(new AlignmentSequence(trainingAlignmentFiles[assignedID], seqNames));
      scrf->setSequences(trainingAlignments.back(), NULL, NULL);
      GTF singleGTF(trainingGTFFiles[assignedID]);
      singleGTF.discardOverlaps();

      //if our CRF does not support 5' or 3' UTRs, delete those annotations from the GTFs
      if (scrf->getStateID("FivePrimeUTR") == -1 && scrf->getStateID("FivePrimeUTRSingle") == -1)
	singleGTF.removeFeatureType("5UTR");
      if (scrf->getStateID("ThreePrimeUTR") == -1)
	singleGTF.removeFeatureType("3UTR");
      
      trainingSegmentations.push_back(Segmentation(singleGTF, *scrf));
      trainingSegmentations.back().setLabels();
    }
    else {
      testingAlignments.push_back(new AlignmentSequence(testingAlignmentFiles[assignedID], seqNames));
      scrf->setSequences(testingAlignments.back(), NULL, NULL);
      GTF singleGTF(testingGTFFiles[assignedID]);
      singleGTF.discardOverlaps(); 

      //if our CRF does not support 5' or 3' UTRs, delete those annotations from the GTFs
      if (scrf->getStateID("FivePrimeUTR") == -1 && scrf->getStateID("FivePrimeUTRSingle") == -1)
	singleGTF.removeFeatureType("5UTR");
      if (scrf->getStateID("ThreePrimeUTR") == -1)
	singleGTF.removeFeatureType("3UTR");

      testingSegmentations.push_back(Segmentation(singleGTF, *scrf));
      testingSegmentations.back().setLabels();
    }
  }
}
