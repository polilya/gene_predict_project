#ifndef SVMTRAINER_H
#define SVMTRAINER_H

#include "SemiCRF.h"

enum svm_command_t {RECEIVE_SVM_TRAINING_DATA,
		    CV_SVM,
		    SEND_SVM,
		    RECEIVE_SVM,
		    SVM_TRAINING_COMPLETE,
		    ABORT_SVM_TRAINING,
		    DO_NOTHING,
		    READY};

class SVMTrainer {
 public:
  int id;
  int numProcs;
  SemiCRF* scrf;

  vector<AlignmentSequence*> trainingAlignments;
  vector<AlignmentSequence*> testingAlignments;
  vector<Segmentation> trainingSegmentations;
  vector<Segmentation> testingSegmentations;
  vector<string> trainingAlignmentFiles;
  vector<string> trainingGTFFiles;
  vector<string> testingAlignmentFiles;
  vector<string> testingGTFFiles;

  SVMTrainer(SemiCRF* scrf, string trainingList, string testingList);
  ~SVMTrainer();
  void train();
  void collectExamples(AlignmentSequence* alignment, Segmentation& segmentation, bool training);
  void addNoncodingDecoy(SVMFS* svmfs, AlignmentSequence* alignment, Segmentation& segmentation, bool training);
  void addCloseDecoy(SVMFS* svmfs, AlignmentSequence* alignment, 
		     int transitionID, pos_t pos, bool training);
  void setWorkAssignment(vector<int> assignment, bool training);

#ifdef MULTI
  void sendCommand(svm_command_t command, int destinationID);
  void broadcastCommand(svm_command_t command);
  void sendSparseVectors(int destinationID, vector<vector<int> >& sparseVectors);
  void receiveSparseVectors(int sourceID, vector<vector<int> >& sparseVectors);
#endif

};

void* cvThread(void* svmfs);

#endif
