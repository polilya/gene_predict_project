#include "SVMFS.h"
#include <stdlib.h>
#include <iomanip>
#include <pthread.h>

SVMFS::SVMFS(int id, svmfs_t svmfsType, kernel_t kernelType, int kernelOrder, int numBins, int length,
	     int offset, weight_t sampleRate, vector<string>& seqNames) {
  svm = NULL;
  this->id = id;
  this->type = SVM;
  this->svmfsType = svmfsType;
  this->kernelType = kernelType;
  this->kernelOrder = kernelOrder;
  this->numBins = numBins;
  this->length = length;
  this->offset = offset;
  this->sampleRate = sampleRate;
  this->seqNames = seqNames;
  cvAccuracy = -1.0;
  cvComplete = false;

  if (svmfsType == CODON) {
    if (length % 3 != 0)
      fatalError("Length of Codon SVMFS is not a multiple of 3");
    setCodonFeatureTable();
  }

  frame = -1;

  int numValues = numBins;
  valueToWeight.resize(numValues);
  valueToGlobalParamID.resize(numValues);
  valueToParam.resize(numValues);
  for (int i=0; i<numBins; i++)
    valueToParam[i] = i;
}

void SVMFS::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;

  seqIDs.clear();
  for (int i=0; i<seqNames.size(); i++) {
    seqIDs.push_back(alignment->getSpeciesID(seqNames[i]));
  }

  if (numericSeq != NULL) {
    bins[PLUS] = bins[NONE] = numericSeq->tags[2 * id];
    bins[MINUS] = numericSeq->tags[2 * id + 1];
    coeffs[PLUS] = coeffs[NONE] = numericSeq->sequenceArray[2 * id];
    coeffs[MINUS] = numericSeq->sequenceArray[2 * id + 1];
  }
}

weight_t SVMFS::calcDecisionValue(vector<int>& x) {

  if (kernelType == POLYNOMIAL) {
    /*
    if (kernelOrder == 2) {
      //use weight vector directly
      int numFeatures = length * seqNames.size() * DNA_CHARS;
      float decision = 0;
      vector<int>::iterator indicesStart = x.indices.begin();
      vector<int>::iterator indicesEnd = x.indices.end();
      vector<int>::iterator firstIndexPtr;
      vector<int>::iterator secondIndexPtr;
      for (firstIndexPtr = indicesStart; firstIndexPtr != indicesEnd; ++firstIndexPtr) {
	int offset = *firstIndexPtr * numFeatures;
	for (secondIndexPtr = firstIndexPtr; secondIndexPtr != indicesEnd; ++secondIndexPtr) {
	  decision += weights[offset + *secondIndexPtr];
	}
      }
      decision -= *(svm->rho);
      return decision;
    }
    else {
    */
      // Use svLookup
      vector<int> dotProducts(svm->l, 0);
      for (int i=0; i<x.size(); i++) {
	int index = x[i];
	//for (int j=0; j<svLookup[index].size(); j++)
	//dotProducts[svLookup[index][j]]++;
	for (vector<int>::iterator iter = svLookup[index].begin(); iter != svLookup[index].end(); ++iter)
	  dotProducts[*iter]++;
      }
      
      weight_t decisionValue = 0;
      for (int i=0; i<svm->l; i++) {
	weight_t contribution = svm->sv_coef[0][i];
	for (int j=0; j<kernelOrder; j++) 
	  contribution *= dotProducts[i];
	decisionValue += contribution;
      }
      decisionValue -= *(svm->rho);
      
      return decisionValue;
      //}
  }
  else {
    // Use libsvm
    struct svm_node* nodes = (struct svm_node*)malloc ((x.size()+1) * sizeof(struct svm_node));
    for (int i=0; i<x.size(); i++) {
      nodes[i].index = x[i];
      nodes[i].value = 1;
    }
    nodes[x.size()].index = -1;  //sentinel node
    double decisionValue;
    svm_predict_values(svm, nodes, &decisionValue);
    
    free(nodes);
    
    return decisionValue;
  }
}

void SVMFS::setWeightVector() {
  if (kernelType != POLYNOMIAL || kernelOrder != 2)
    fatalError("setWeightVector only supports quadratic kernels currently");

  int numFeatures = length * seqNames.size() * DNA_CHARS;
  weights.clear();
  weights.resize(numFeatures * numFeatures, 0.0);
  
  for (int i=0; i<svm->l; i++) {
    for (int j=0; svm->SV[i][j].index != -1; j++) {
      weights[numFeatures*svm->SV[i][j].index + svm->SV[i][j].index] += svm->sv_coef[0][i];
      for (int k=j+1; svm->SV[i][k].index != -1; k++) {
	//count both (j,k) and (k,j) in a single cell to save time later
	weights[numFeatures*svm->SV[i][j].index + svm->SV[i][k].index] += 2 * svm->sv_coef[0][i];
      }
    }
  }
}

void SVMFS::setCodonFeatureTable() {
  numCodonFeatures = 75;
  int numTrimers = DNA_CHARS * DNA_CHARS * DNA_CHARS;

  codonFeatures.resize(numTrimers);
  for (int charOneIndex=0; charOneIndex<DNA_CHARS; charOneIndex++) { 
    for (int charTwoIndex=0; charTwoIndex<DNA_CHARS; charTwoIndex++) { 
      for (int charThreeIndex=0; charThreeIndex<DNA_CHARS; charThreeIndex++) {

	char charOne = INDEX_TO_DNA[charOneIndex];
	char charTwo = INDEX_TO_DNA[charTwoIndex];
	char charThree = INDEX_TO_DNA[charThreeIndex];
 
	int trimerIndex = DNA_CHARS * DNA_CHARS * charOneIndex + 
	  DNA_CHARS * charTwoIndex + charThreeIndex;
	int featureIndex;
	
	if (charOne == 'N' || charTwo == 'N' || charThree == 'N') {
	  //special case, codon contains unusual character
	  featureIndex = 75;
	}
	else if (charOne == '.' || charTwo == '.' || charThree == '.') {
	  //codon contains an unaligned position
	  if (charOne == '.' && charTwo == '.' && charThree == '.')
	    featureIndex = 74;  //whole codon is unaligned
	  else if (charOne == '.' && charTwoIndex < 4 && charThreeIndex < 5)
	    featureIndex = 72;  //alignment start
	  else if (charOne == '.' && charTwo == '.' && charThreeIndex < 4)
	    featureIndex = 72;  //alignment start
	  else if (charThree == '.' && charTwoIndex < 4 && charOneIndex < 5)
	    featureIndex = 73;  //alignment end
	  else if (charThree == '.' && charTwo == '.' && charOneIndex < 4)
	    featureIndex = 73;  //alignment end
	  else
	    featureIndex = 75;  //something weird with this codon
	}
	else if (charOne == '_' || charTwo == '_' || charThree == '_') {
	  //codon contains gap(s) but no unaligned positions
	  int gapIndex = 0;
	  if(charOne == '_')
	    gapIndex += 4;
	  if (charTwo == '_')
	    gapIndex += 2;
	  if (charThree == '_')
	    gapIndex += 1;
	  featureIndex = 64 + gapIndex;
	}
	else {
	  //codon contains only base pairs
	  featureIndex = 16 * charOneIndex + 4 * charTwoIndex + charThreeIndex + 1;
	}
 
	codonFeatures[trimerIndex] = featureIndex;

	//cerr << trimerIndex << "\t" << charOne << charTwo << charThree << "\t" << featureIndex << endl;
      }
    }
  }
}

bool SVMFS::buildInputVector(vector<int>& x, pos_t pos, strand_t strand) {
  x.clear();

  if (svmfsType == POSITIONAL) {
    pos_t start;
    pos_t row;
    
    if (strand == PLUS || strand == NONE)
      start = pos - offset;
    else
      start = alignment->length - 1 - pos - offset + 1;  //extra +1 for transitions here!
    
    if (start < 0 || start + length > alignment->length)
      return false;  //out of bounds, don't build a vector
    
    int offset = 0;
    
    for (int i=0; i<seqIDs.size(); i++) {
      int row;
      if (strand == PLUS || strand == NONE)
	row = 2 * seqIDs[i];
      else
	row = 2 * seqIDs[i] + 1;
      
      for (pos_t j=start; j<start + length; j++) {
	x.push_back(offset + DNA_TO_INDEX[alignment->sequenceArray[row][j]]);
	offset += DNA_CHARS;
      }
    }
    return true;  //success
  }
  else if (svmfsType == CODON) {
    pos_t start;
    pos_t row;
    
    if (strand == PLUS || strand == NONE)
      start = pos - offset;
    else
      start = alignment->length - 1 - pos - offset;
    
    if (start < 0 || start + length > alignment->length)
      return false;  //out of bounds, don't build a vector
    
    int offset = 0;
    
    for (int i=0; i<seqIDs.size(); i++) {
      int row;
      if (strand == PLUS || strand == NONE)
	row = 2 * seqIDs[i];
      else
	row = 2 * seqIDs[i] + 1;
      
      for (pos_t j=start; j<start + length; j += 3) {
	int charOneIndex = DNA_TO_INDEX[alignment->sequenceArray[row][j]];
	int charTwoIndex = DNA_TO_INDEX[alignment->sequenceArray[row][j+1]];
	int charThreeIndex = DNA_TO_INDEX[alignment->sequenceArray[row][j+2]];
	int trimerIndex = DNA_CHARS * DNA_CHARS * charOneIndex + 
	  DNA_CHARS * charTwoIndex + charThreeIndex;
	x.push_back(offset + codonFeatures[trimerIndex]);
	offset += numCodonFeatures;
      }
    }
    return true;  //success    
  }

  return false;
}

void SVMFS::buildProblem(struct svm_problem& problem, bool training) {
  vector<vector<int> >* trueExamples;
  vector<vector<int> >* decoys;
  if (training) {
    trueExamples = &trainingTrue;
    decoys = &trainingDecoys;
  }
  else {
    trueExamples = &testingTrue;
    decoys = &testingDecoys;
  }

  problem.l = trueExamples->size() + decoys->size();

  problem.y = (double*) malloc (problem.l * sizeof(double));
  problem.x = (struct svm_node**) malloc (problem.l * sizeof(struct svm_node*));
  int probIndex = 0;
  for (int i=0; i<trueExamples->size(); i++) {
    problem.x[probIndex] = (struct svm_node*) malloc (((*trueExamples)[i].size()+1) * sizeof(struct svm_node));
    for (int j=0; j<(*trueExamples)[i].size(); j++) {
      problem.x[probIndex][j].index = (*trueExamples)[i][j];
      problem.x[probIndex][j].value = 1;
    }
    problem.x[probIndex][(*trueExamples)[i].size()].index = -1;  //sentinel node
    problem.y[probIndex] = 1;
    probIndex++;
  }
  for (int i=0; i<decoys->size(); i++) {
    problem.x[probIndex] = (struct svm_node*) malloc (((*decoys)[i].size()+1) * sizeof(struct svm_node));
    for (int j=0; j<(*decoys)[i].size(); j++) {
      problem.x[probIndex][j].index = (*decoys)[i][j];
      problem.x[probIndex][j].value = 1;
    }
    problem.x[probIndex][(*decoys)[i].size()].index = -1;  //sentinel node
    problem.y[probIndex] = -1;
    probIndex++;
  }
}

void SVMFS::freeProblem(struct svm_problem& problem) {
  for (int i=0; i<problem.l; i++) {
    free(problem.x[i]);
  }

  free(problem.x);
  free(problem.y);
}

void SVMFS::setSVMTrainingParameters(struct svm_parameter& parameters, weight_t C, bool probabilities) {
  parameters.svm_type = C_SVC;
  if (kernelType == POLYNOMIAL) {
    parameters.kernel_type = POLY;
    parameters.degree = kernelOrder;
  }
  else if (kernelType == GAUSSIAN) {
    parameters.kernel_type = RBF;
    parameters.degree = 3;
  }
  parameters.gamma = 1;
  parameters.coef0 = 0;
  parameters.nu = 0.5;
  parameters.cache_size = 1000;
  parameters.C = C;
  parameters.eps = 1e-3;
  parameters.p = 0.1;
  parameters.shrinking = 1;
  
  if (probabilities)
    parameters.probability = 1;
  else
    parameters.probability = 0;

  parameters.nr_weight = 0;
  parameters.weight_label = NULL;
  parameters.weight = NULL;
}

svm_t* SVMFS::copySVM(svm_t* svm) {
  svm_t* copy = (struct svm_model*) malloc (sizeof(struct svm_model));
  copy->l = svm->l;
  copy->label = (int*) malloc (2 * sizeof(int));
  copy->label[0] = svm->label[0];
  copy->label[1] = svm->label[1];
  copy->nSV = (int*) malloc (2 * sizeof(int));
  copy->nSV[0] = svm->nSV[0];
  copy->nSV[1] = svm->nSV[1];
  copy->free_sv = 1;
  copy->nr_class = 2;
  copy->param.svm_type = C_SVC;
  if (kernelType == POLYNOMIAL) {
    copy->param.kernel_type = POLY;
    copy->param.degree = kernelOrder;
  }
  else if (kernelType == GAUSSIAN)
    copy->param.kernel_type = GAUSSIAN;
  copy->param.gamma = 1;
  copy->param.coef0 = 0;
  copy->rho = (weight_t*) malloc (sizeof(weight_t));
  *(copy->rho) = *(svm->rho);
  copy->probA = (weight_t*) malloc (sizeof(weight_t));
  *(copy->probA) = 0;
  copy->probB = (weight_t*) malloc (sizeof(weight_t));
  *(copy->probB) = 0;

  copy->sv_coef = (weight_t**) malloc (sizeof(weight_t*));
  copy->sv_coef[0] = (weight_t*) malloc (svm->l * sizeof(weight_t));
  for (int i=0; i<svm->l; i++)
    copy->sv_coef[0][i] = svm->sv_coef[0][i];

  copy->SV = (struct svm_node**) malloc (svm->l * sizeof(struct svm_node*));
  for (int i=0; i<svm->l; i++) {
    int entries = 0;
    for (int j=0; svm->SV[i][j].index != -1; j++)
      entries++;
    copy->SV[i] = (struct svm_node*) malloc ((entries + 1) * sizeof(struct svm_node));
    for (int j=0; j<=entries; j++)
      copy->SV[i][j] = svm->SV[i][j];
  }

  return copy;
}

weight_t SVMFS::cvSVM() {
  //assume C has already been set!

  //if (svm != NULL)
  //svm_destroy_model(svm);

  struct svm_problem problem;
  struct svm_parameter parameters;
  
  cvComplete = false;
  cvAccuracy = -1.0;

  //build problem and train SVM
  buildProblem(problem, true);
  setSVMTrainingParameters(parameters, C, false);
  svm_t* tempSVM = svm_train(&problem, &parameters);
  svm = copySVM(tempSVM);
  svm_destroy_model(tempSVM);
  freeProblem(problem);

  //set lookup tables for fast classification
  setSVLookup();
  //if (kernelType == POLYNOMIAL && kernelOrder == 2)
  //setWeightVector();

  //set bin limits
  vector<weight_t> decisionValues;
  for (int i=0; i<testingTrue.size(); i++)
    decisionValues.push_back(calcDecisionValue(testingTrue[i]));
  sort(decisionValues.begin(), decisionValues.end());
  weight_t lowLimit = decisionValues[(int)(0.01 * decisionValues.size())];
  weight_t highLimit = decisionValues[(int)(0.99 * decisionValues.size())];
  binLimits.resize(numBins);
  binLimits[0] = LOG_ZERO;  //negative infinity
  for (int i=1; i<numBins; i++)
    binLimits[i] = lowLimit + (i-1) * (highLimit - lowLimit)/(numBins-2.0);

  //get generalization accuracy from testing data and record completion
  cvAccuracy = getAccuracy(false);
  cvComplete = true;

  return cvAccuracy;
}

weight_t SVMFS::getAccuracy(bool training) {

  struct svm_problem problem;
  buildProblem(problem, training);
  
  int correct = 0; 
  for (int i=0; i<problem.l; i++) {
    if (svm_predict(svm, problem.x[i]) == problem.y[i])
      correct++;
  }
  
  weight_t accuracy = (weight_t)correct / (weight_t)problem.l;

  freeProblem(problem);

  return accuracy;
}

void SVMFS::testCalcDecisionValue() {
  struct svm_problem problem;
  buildProblem(problem, false);

  for (int i=0; i<problem.l; i++) {
    double libsvmDecisionValue;
    weight_t decisionValue;
    svm_predict_values(svm, problem.x[i], &libsvmDecisionValue);
    if (i < testingTrue.size())
      decisionValue = calcDecisionValue(testingTrue[i]);
    else
      decisionValue = calcDecisionValue(testingDecoys[i - testingTrue.size()]);
    cerr << "libsvm: " << libsvmDecisionValue << ", optimized: " << decisionValue << endl;
  }

  freeProblem(problem);
}

//set initial guess based on the testing data
void SVMFS::makeInitialGuess() {

  vector<weight_t> decoyBins(numBins, 0);  //records where decoy decision values fall
  vector<weight_t> trueBins(numBins, 0);   //records where true decision values fall
  weight_t trueDenominator = 0;
  weight_t decoyDenominator = 0;

  for (int i=0; i<testingTrue.size(); i++) {
    weight_t decisionValue = calcDecisionValue(testingTrue[i]);
    int bin = lower_bound(binLimits.begin(), binLimits.end(), decisionValue) - binLimits.begin() - 1;
    if (bin != 0 && bin != binLimits.size() - 1) {
      weight_t coeff = (binLimits[bin+1] - decisionValue) / (binLimits[bin+1] - binLimits[bin]);
      trueBins[bin] += coeff;
      trueBins[bin+1] += 1.0 - coeff;
      trueDenominator += 1.0;
    }
  }
  for (int i=0; i<testingDecoys.size(); i++) {
    weight_t decisionValue = calcDecisionValue(testingDecoys[i]);
    int bin = lower_bound(binLimits.begin(), binLimits.end(), decisionValue) - binLimits.begin() - 1;
    if (bin != 0 && bin != binLimits.size() - 1) {
      weight_t coeff = (binLimits[bin+1] - decisionValue) / (binLimits[bin+1] - binLimits[bin]);
      decoyBins[bin] += coeff;
      decoyBins[bin+1] += 1.0 - coeff;
      decoyDenominator += 1.0;
    }
  }

  //add +1 pseudocounts
  for (int i=0; i<numBins; i++) {
    trueBins[i] += 1.0;
    decoyBins[i] += 1.0;
    trueDenominator += 1.0;
    decoyDenominator += 1.0;    
  }

  cerr << "Initial guess for SVM id " << id << endl;
  for (int i=0; i<numBins; i++) {  //calculate log probability ratio
    cerr << "Bin " << i << ": " << trueBins[i] << " true examples of " << trueDenominator << " fell into this bin..." 
	 << decoyBins[i] << " decoy examples of " << decoyDenominator << " fell into this bin" << endl; 
    weight_t logRatio = log( (trueBins[i] / trueDenominator) / (decoyBins[i] / decoyDenominator) );
    paramToWeight[i] = logRatio;
  }
}

//add a random decoy from the alignment matching one of the k-mer
//requirements of the given transition
void SVMFS::addRandomDecoy(TransitionFeature& transition, bool training) {
  //randomly select a required kmer
  char* requiredKmer;
  int adjustment = 0;
  if (transition.toStart.size() == 0) {
    int requiredIndex = (int)((double)rand() / RAND_MAX * transition.fromEnd.size());
    requiredKmer = transition.fromEnd[requiredIndex];
    adjustment = strlen(requiredKmer);
  }
  else if (transition.fromEnd.size() == 0) {
    int requiredIndex = (int)((double)rand() / RAND_MAX * transition.toStart.size());
    requiredKmer = transition.toStart[requiredIndex];
  }
  else {
    fatalError("SVMFS: required k-mer array sizes incorrect");
  }

  int kmerLength = strlen(requiredKmer);
  vector<int> x;
  pos_t j;
  while (1) {
    j = (int)((weight_t)rand() / (weight_t)RAND_MAX * (weight_t)alignment->length);
    if (strncmp(alignment->sequenceArray[0] + j, requiredKmer, kmerLength) == 0 &&
	buildInputVector(x, j + adjustment, transition.strand))
      break;
  }

  if (training)
    trainingDecoys.push_back(x);
  else
    testingDecoys.push_back(x);
}


void SVMFS::readSVM(ifstream& svmStream) {
  string s1 = "";
  string s2;
  string s3;

  svm = (struct svm_model*) malloc (sizeof(struct svm_model));

  //read SVM parameters
  while (s1 != "SV") {
    svmStream >> s1;

    if (s1 == "svm_type") {
      svmStream >> s2;
      if (s2 == "c_svc")
	svm->param.svm_type = C_SVC;
      else
	fatalError("Unknown SVM type: " + s2);
    }
    else if (s1 == "kernel_type") {
      svmStream >> s2;
      if (s2 == "polynomial") {
	svm->param.kernel_type = POLY;
      }
      else
	fatalError("Unknown kernel type: " + s2);
    }
    else if (s1 == "degree") {
      svmStream >> svm->param.degree;
    }
    else if (s1 == "gamma") {
      svmStream >> svm->param.gamma;
    }
    else if (s1 == "coef0") {
      svmStream >> svm->param.coef0;
    }
    else if (s1 == "nr_class") {
      svmStream >> svm->nr_class;
      svm->label = (int*) malloc (svm->nr_class * sizeof(int));
      svm->nSV = (int*) malloc (svm->nr_class * sizeof(int));
    }
    else if (s1 == "total_sv") {
      svmStream >> svm->l;

      svm->sv_coef = (double**) malloc (sizeof(double*));
      svm->sv_coef[0] = (double*) malloc (svm->l * sizeof(double));
      svm->SV = (struct svm_node**) malloc (svm->l * sizeof(struct svm_node*));
    }
    else if (s1 == "rho") {
      svm->rho = (double*) malloc (svm->nr_class * (svm->nr_class-1) / 2 * sizeof(double));
      for (int i=0; i<svm->nr_class * (svm->nr_class-1) / 2; i++)
	svmStream >> svm->rho[i];
    }
    else if (s1 == "label") {
      for (int i=0; i<svm->nr_class; i++)
	svmStream >> svm->label[i];
    }
    else if (s1 == "nr_sv") {
      for (int i=0; i<svm->nr_class; i++)
	svmStream >> svm->nSV[i];
    }
    else if (s1 == "SV") {
      //we're done here

    }
    else
      fatalError("Unknown SVM parameter tag: " + s1);
  }

  //read support vectors and their weights
  for (int i=0; i<svm->l; i++) {
    svmStream >> svm->sv_coef[0][i];
    char svLine[16384];
    svmStream.getline(svLine, 16384);
    string svLineStr(svLine);
    stringstream ss(svLineStr);
    vector<string> entries;
    string entry;
    while (ss >> entry)
      entries.push_back(entry);
    svm->SV[i] = (struct svm_node*) malloc ((entries.size() + 1) * sizeof(struct svm_node));
    for (int j=0; j<entries.size(); j++) {
      int colonPos = entries[j].find(':', 0);
      svm->SV[i][j].index = atoi(entries[j].substr(0, colonPos).c_str());
      svm->SV[i][j].value = atof(entries[j].substr(colonPos+1, entries[j].length() - colonPos - 1).c_str());
    }
    svm->SV[i][entries.size()].index = -1;  //sentinel value
  }

  setSVLookup();
  //if (kernelType == POLYNOMIAL && kernelOrder == 2)
  //setWeightVector();
}

void SVMFS::printSVM(ostream& svmStream) {
  if (svm == NULL) {
    svmStream << "No SVM" << endl;
    return;
  }

  const svm_parameter& param = svm->param;

  svmStream << "svm_type " << svm_type_table[param.svm_type] << endl;
  svmStream << "kernel_type " << kernel_type_table[param.kernel_type] << endl;
  
  if(param.kernel_type == POLY)
    svmStream << "degree " << param.degree << endl;
  
  if(param.kernel_type == POLY || param.kernel_type == RBF || param.kernel_type == SIGMOID)
    svmStream << "gamma " << param.gamma << endl;
  
  if(param.kernel_type == POLY || param.kernel_type == SIGMOID)
    svmStream << "coef0 " << param.coef0 << endl;
  
  int nr_class = svm->nr_class;
  int l = svm->l;
  svmStream << "nr_class " << nr_class << endl;
  svmStream << "total_sv " << l << endl;
  
  {
    svmStream << "rho";
    for(int i=0;i<nr_class*(nr_class-1)/2;i++)
      svmStream << " " << svm->rho[i];
    svmStream << endl;
  }
  
  if(svm->label)
    {
      svmStream << "label";
      for(int i=0;i<nr_class;i++)
	svmStream << " " << svm->label[i];
      svmStream << endl;
    }
  
  if(svm->nSV)
    {
      svmStream << "nr_sv";
      for(int i=0;i<nr_class;i++)
	svmStream << " " << svm->nSV[i];
      svmStream << endl;
    }

    svmStream << "SV" << endl;
    const double * const *sv_coef = svm->sv_coef;
    const svm_node * const *SV = svm->SV;
  
    for(int i=0;i<l;i++)
      {
	for(int j=0;j<nr_class-1;j++)
	  svmStream << sv_coef[j][i] << " ";
	
	const svm_node *p = SV[i];
	while(p->index != -1)
	  {
	    svmStream << p->index << ":" << p->value << " ";
	    p++;
	  }
	svmStream << endl;
      }
}

void SVMFS::printAlignmentBlock(ostream& os, vector<int>& x) {
  int offset = 0;
  for (int i=0; i<x.size(); i++) {
    os << INDEX_TO_DNA[x[i] - offset];
    offset += DNA_CHARS;
    if ((i+1) % length == 0)
      os << endl;
  }
}

void SVMFS::setSVLookup() {
  if (svmfsType == POSITIONAL) {
    svLookup.resize(length * seqNames.size() * DNA_CHARS);
  }
  else if (svmfsType == CODON) {
    svLookup.resize(length/3 * seqNames.size() * codonFeatures.size());
  }
  for (int i=0; i<svLookup.size(); i++)
    svLookup[i].clear();

  for (int i=0; i<svm->l; i++) {
    for (int j=0; svm->SV[i][j].index != -1; j++) {
      svLookup[svm->SV[i][j].index].push_back(i);
    }
  }
}

#ifdef MULTI

void SVMFS::sendSVM(int destinationID) {

  //calculate total number of sparse vector entries in SVM model
  int numValues = 0;
  for (int i=0; i<svm->l; i++) {
    for (int j=0; svm->SV[i][j].index != -1; j++)
      numValues++;
  }

  int bufferOne[7];
  bufferOne[0] = svm->l;
  bufferOne[1] = numValues;
  bufferOne[2] = svm->label[0];
  bufferOne[3] = svm->label[1];
  bufferOne[4] = svm->nSV[0];
  bufferOne[5] = svm->nSV[1];
  bufferOne[6] = svm->free_sv;
  MPI_Send(bufferOne, 7, MPI_INT, destinationID, 0, MPI_COMM_WORLD);

  weight_t* bufferTwo = (weight_t *) malloc ((3 + numBins + svm->l) * sizeof(weight_t));
  bufferTwo[0] = *(svm->rho);
  bufferTwo[1] = 0;
  bufferTwo[2] = 0;
  //bufferTwo[1] = *(svm->probA);
  //bufferTwo[2] = *(svm->probB);
  for (int i=0; i<numBins; i++)
    bufferTwo[3 + i] = binLimits[i];
  for (int i=0; i<svm->l; i++)
    bufferTwo[3 + numBins + i] = svm->sv_coef[0][i];
  MPI_Send(bufferTwo, 3 + numBins + svm->l, MPI_WEIGHT_T, destinationID, 0, MPI_COMM_WORLD);

  int* bufferThree = (int *) malloc ((numValues + svm->l) * sizeof(int));
  int bufferIndex = 0;
  for (int i=0; i<svm->l; i++) {
    for (int j=0; svm->SV[i][j].index != -1; j++)
      bufferThree[bufferIndex++] = svm->SV[i][j].index;
    bufferThree[bufferIndex++] = -1;
  }
  MPI_Send(bufferThree, numValues + svm->l, MPI_INT, destinationID, 0, MPI_COMM_WORLD);

  weight_t* bufferFour = (weight_t *) malloc (numValues * sizeof(weight_t));
  bufferIndex = 0;
  for (int i=0; i<svm->l; i++) {
    for (int j=0; svm->SV[i][j].index != -1; j++)
      bufferFour[bufferIndex++] = svm->SV[i][j].value;
  }
  MPI_Send(bufferFour, numValues, MPI_WEIGHT_T, destinationID, 0, MPI_COMM_WORLD);

  free(bufferTwo);
  free(bufferThree);
  free(bufferFour);
}

void SVMFS::receiveSVM(int sourceID) {

  if (svm != NULL)
    svm_destroy_model(svm);
  svm = (struct svm_model*) malloc (sizeof(struct svm_model));
    
  int bufferOne[7];
  MPI_Recv(bufferOne, 7, MPI_INT, sourceID, 0, MPI_COMM_WORLD, NULL);
  svm->l = bufferOne[0];
  int numValues = bufferOne[1];
  svm->label = (int*) malloc (2 * sizeof(int));
  svm->label[0] = bufferOne[2];
  svm->label[1] = bufferOne[3];
  svm->nSV = (int*) malloc (2 * sizeof(int));
  svm->nSV[0] = bufferOne[4];
  svm->nSV[1] = bufferOne[5];
  svm->free_sv = bufferOne[6];

  svm->nr_class = 2;
  svm->param.svm_type = C_SVC;
  if (kernelType == POLYNOMIAL) {
    svm->param.kernel_type = POLY;
    svm->param.degree = kernelOrder;
  }
  else if (kernelType == GAUSSIAN)
    svm->param.kernel_type = GAUSSIAN;
  svm->param.gamma = 1;
  svm->param.coef0 = 0;

  weight_t* bufferTwo = (weight_t *) malloc ((3 + numBins + svm->l) * sizeof(weight_t));
  MPI_Recv(bufferTwo, 3 + numBins + svm->l, MPI_WEIGHT_T, sourceID, 0, MPI_COMM_WORLD, NULL);
  svm->rho = (weight_t*) malloc (sizeof(weight_t));
  *(svm->rho) = bufferTwo[0];
  svm->probA = (weight_t*) malloc (sizeof(weight_t));
  *(svm->probA) = bufferTwo[1];
  svm->probB = (weight_t*) malloc (sizeof(weight_t));
  *(svm->probB) = bufferTwo[2];
  binLimits.resize(numBins);
  for (int i=0; i<numBins; i++)
    binLimits[i] = bufferTwo[3 + i];
  svm->sv_coef = (weight_t**) malloc (sizeof(weight_t*));
  svm->sv_coef[0] = (weight_t*) malloc (svm->l * sizeof(weight_t));
  for (int i=0; i<svm->l; i++)
    svm->sv_coef[0][i] = bufferTwo[3 + numBins + i];
  

  int* bufferThree = (int *) malloc ((numValues + svm->l) * sizeof(int));
  MPI_Recv(bufferThree, numValues + svm->l, MPI_INT, sourceID, 0, MPI_COMM_WORLD, NULL);  
  int bufferIndex = 0;
  svm->SV = (struct svm_node**) malloc (svm->l * sizeof(struct svm_node*));
  for (int i=0; i<svm->l; i++) {
    int entries = 0;
    for (int j=0; bufferThree[bufferIndex + j] != -1; j++)
      entries++;
    svm->SV[i] = (struct svm_node*) malloc ((entries + 1) * sizeof(struct svm_node));
    for (int j=0; bufferThree[bufferIndex + j] != -1; j++)
      svm->SV[i][j].index = bufferThree[bufferIndex + j];
    svm->SV[i][entries].index = -1;
    bufferIndex += entries + 1;
  }
  

  weight_t* bufferFour = (weight_t *) malloc (numValues * sizeof(weight_t));
  MPI_Recv(bufferFour, numValues, MPI_WEIGHT_T, sourceID, 0, MPI_COMM_WORLD, NULL);
  bufferIndex = 0;
  for (int i=0; i<svm->l; i++) {
    for (int j=0; svm->SV[i][j].index != -1; j++)
      svm->SV[i][j].value = bufferFour[bufferIndex++];
  }

  free(bufferTwo);
  free(bufferThree);
  free(bufferFour);

  setSVLookup();

  //if (kernelType == POLYNOMIAL && kernelOrder == 2)
  //setWeightVector();
} 

#endif
