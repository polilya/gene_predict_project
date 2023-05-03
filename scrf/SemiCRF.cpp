#include "SemiCRF.h"
#include <sstream>
#include <map>
#include <set>

#include "MemoryDebug.h"

SemiCRF::SemiCRF(string parameterFile) {
  alignment = NULL;
  numericSeq = NULL;
  estSeq = NULL;

  /* read the parameters */
  readParameters(parameterFile);
  setParameterMap();
  //cerr << "Model has " << parameterMap.size() << " total parameters" << endl;

  /* precompute some things to speed up computation */
  setAllowedLookup();
  setAllShortcuts();

  /* initialize properties of the dynamic programming matrix */
  dpMatrix.numberOfStates = states.size();
  dpMatrix.numberOfImplicitStates = numberOfImplicitStates;
  dpMatrix.numberOfTransitions = transitions.size();
  dpMatrix.numberOfTransitionTypes = transitionTypes.size();
}

SemiCRF::~SemiCRF() {
  for (int i=0; i<sequenceFeatureSets.size(); i++)
    delete sequenceFeatureSets[i];

  freeBoolMatrix(allowedLookup, intpow(DNA_CHARS,allowedKmerLength), transitions.size());
}

stateid_t SemiCRF::getStateID(string s){
  for (int i = 0; i < states.size(); i++) {
    if (states[i].name == s) {
      return i;  
    }
  }
  return -1;
}

stateid_t SemiCRF::getStateTypeID(string s){
  for (int i = 0; i < stateTypes.size(); i++)
    if (stateTypes[i].name == s) return i;  
  return -1;
}

int SemiCRF::getTransitionTypeID(string s) {
  for (int i=0; i<transitionTypes.size(); i++)
    if (transitionTypes[i].name == s) return i;
  return -1;
}

SVMFS* SemiCRF::getSVMFSByID(int id) {
  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    SequenceFeatureSet* sfs = sequenceFeatureSets[i];
    if (sfs->type == SVM) {
      SVMFS* svmfs = (SVMFS*)sfs;
      if (svmfs->id == id)
	  return svmfs;
    }
  }

  return NULL;
}

void SemiCRF::readParameters(string filename) {
  ifstream parameterStream(filename.c_str());

  if (parameterStream.bad())
    fatalError("Could not open parameter file");
    
  readStates(parameterStream);
  readTransitions(parameterStream);
  readLengths(parameterStream);
  readMaskingFeatures(parameterStream);
  readDNAKmerFeatures(parameterStream);
  readDNAKmerPairFeatures(parameterStream);
  readDNAKmerArrayFeatures(parameterStream);
  readDNAKmerPairArrayFeatures(parameterStream);
  readSVMFeatures(parameterStream);
  readESTPositionFeatures(parameterStream);
  readESTTransitionFeatures(parameterStream);
}

void SemiCRF::readStates(ifstream& parameterStream) {
  string stateName;
  string stateTypeName;
  string s;

  //loook for the correct header
  parameterStream >> s;
  if (s != "[States]")
    fatalError("Could not find [States] header in parameter file");

  //read in the states
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; put it back and quit */
      /*
      for (int i=0; i<s.length(); i++)
	parameterStream.unget();
      */
      break;
    }

    State state;

    state.name = s;
    parameterStream >> stateTypeName;
    parameterStream >> s;
    if (s == "Yes")
      state.optimize = true;
    else if (s == "No")
      state.optimize = false;
    else
      fatalError("State optimize field must be \"Yes\" or \"No\"");
    parameterStream >> s;
    if (s == "+")
      state.strand = PLUS;
    else if (s == "-")
      state.strand = MINUS;
    else if (s == ".")
      state.strand = NONE;
    else
      fatalError("Illegal strand in parameter file [States] section");

    parameterStream >> state.startWeight;
    
    parameterStream >> s;    
    if (s != ".") {
      istringstream disallowedStream(s);
      while (! disallowedStream.eof()) {
	char disallowedChars[STRLEN];
	disallowedStream.getline(disallowedChars, STRLEN, ',');
	state.disallowedKmers.push_back(strdup(disallowedChars));
      }
    }

    state.typeID = getStateTypeID (stateTypeName);
    
    if (state.typeID == -1){ /* we need to create a new StateType */
      StateType stateType;
      stateType.name = stateTypeName;
      state.typeID = stateTypes.size();
      stateTypes.push_back (stateType);
    }

    states.push_back (state);
  }

  numberOfImplicitStates = states.size();
  numberOfImplicitTypes = stateTypes.size();
}

void SemiCRF::readTransitions(ifstream& parameterStream) {
  string fromState, toState;
  string transitionTypeName;
  string fromEnd, toStart;
  string disallowedFromEnd, disallowedToStart;
  string s;

  //loook for the correct header
  /*
  parameterStream >> s;
  if (s != "[Transitions]")
    fatalError("Could not find [Transitions] header in parameter file");
  */

  //read in the transitions
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; put it back and quit */
      /*
      for (int i=0; i<s.length(); i++)
	parameterStream.unget();
      */
      break;
    }

    TransitionFeature transition;

    fromState = s;
    parameterStream >> toState;
    transition.fromState = getStateID (fromState);
    transition.toState = getStateID (toState);
    if (transition.fromState == -1) {
      cerr << "Illegal transition from state " << fromState << endl;
      fatalError("Illegal transition in parameter file");
    }
    if (transition.toState == -1) {
      cerr << "Illegal transition to state " << toState << endl;
      fatalError("Illegal transition in parameter file");
    }
    if (states[transition.fromState].strand != NONE)
      transition.strand = states[transition.fromState].strand;
    else
      transition.strand = states[transition.toState].strand;

    parameterStream >> transitionTypeName;
    parameterStream >> transition.weight >> fromEnd >> toStart >> disallowedFromEnd >> disallowedToStart;
    if (fromEnd != ".") {
      istringstream fromEndStream(fromEnd);
      while (! fromEndStream.eof()) {
	char requiredChars[STRLEN];
	fromEndStream.getline(requiredChars, STRLEN, ',');
	transition.fromEnd.push_back(strdup(requiredChars));
      }
    }
    if (toStart != ".") {
      istringstream toStartStream(toStart);
      while (! toStartStream.eof()) {
     	char requiredChars[STRLEN];
	toStartStream.getline(requiredChars, STRLEN, ',');
	transition.toStart.push_back(strdup(requiredChars));
      }
    }
    if (disallowedFromEnd != ".") {
      istringstream disallowedFromEndStream(disallowedFromEnd);
      while (! disallowedFromEndStream.eof()) {
     	char disallowedChars[STRLEN];
	disallowedFromEndStream.getline(disallowedChars, STRLEN, ',');
	transition.disallowedFromEnd.push_back(strdup(disallowedChars));
      }
    }
    if (disallowedToStart != ".") {
      istringstream disallowedToStartStream(disallowedToStart);
      while (! disallowedToStartStream.eof()) {
     	char disallowedChars[STRLEN];
	disallowedToStartStream.getline(disallowedChars, STRLEN, ',');
	transition.disallowedToStart.push_back(strdup(disallowedChars));
      }
    }

    transition.typeID = getTransitionTypeID (transitionTypeName);
    
    if (transition.typeID == -1 && transitionTypeName != ".") {  /* we need to create a new TransitionType */
      TransitionType transitionType;
      transitionType.name = transitionTypeName;
      transition.typeID = transitionTypes.size();
      transitionTypes.push_back (transitionType);
    }

    transitions.push_back (transition);
  }
    
  /* tell the states about each transition involving them */
  for (int i=0; i<transitions.size(); i++) {
    states[transitions[i].toState].transitionsTo.push_back(i);
    states[transitions[i].fromState].transitionsFrom.push_back(i);
  }
}


void SemiCRF::readLengths(ifstream& parameterStream) {
  string s;
  int stateTypeID;
  int maxLength;
  string parameterization;

  //loook for the correct header
  /*
  parameterStream >> s;
  if (s != "[Lengths]")
    fatalError("Could not find [Lengths] header in parameter file");
  */  

  //read in the lengths
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; put it back and quit */
      /*
      for (int i=0; i<s.length(); i++)
	parameterStream.unget();
      */
      break;
    }
    
    if (getStateTypeID(s) != -1) {
      stateTypeID = getStateTypeID(s);
      parameterStream >> maxLength;
      parameterStream >> parameterization;
      stateTypes[stateTypeID].lengthFeatures = new LengthFS(maxLength, parameterization);
    }
    else {
      //set this parameter
      stateTypes[stateTypeID].lengthFeatures->paramToWeight.push_back(atof(s.c_str()));
    }
  }
}

void SemiCRF::readMaskingFeatures(ifstream& parameterStream) {
  string s;
  SequenceFeatureSet* currentSFS;
  float param;

  //read in the features
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; quit */
      break;
    }
    
    if (sscanf(s.c_str(), "%f", &param)) {
      //set this parameter
      currentSFS->paramToWeight.push_back(param);
    }
    else { /* read header */ 
      istringstream headerStream(s);
      vector<int> stateTypeIDs;
      while (! headerStream.eof()) {
	char stateTypeName[STRLEN];
	headerStream.getline(stateTypeName, STRLEN, ',');
	int typeID = getStateTypeID(string(stateTypeName));
	if (typeID == -1) {
	  cerr << "Unrecognized state type name: " << stateTypeName << endl;
	  fatalError("Could not parse MaskingFeatures");
	}
        stateTypeIDs.push_back(typeID);
      }
      sequenceFeatureSets.push_back( new MaskingFeature() );
      for (int i=0; i<stateTypeIDs.size(); i++) {
	sequenceFeatureSets.back()->stateTypeNames.push_back(stateTypes[stateTypeIDs[i]].name);
	stateTypes[stateTypeIDs[i]].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
      }
      currentSFS = sequenceFeatureSets.back();
    }
  }
}

void SemiCRF::readDNAKmerFeatures(ifstream& parameterStream) {
  string s;
  int k;
  string frameString;
  int frame;
  string speciesName;
  string parameterization;
  SequenceFeatureSet* currentSFS;
  float param;

  //read in the features
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; quit */
      break;
    }
    
    if (sscanf(s.c_str(), "%f", &param)) {
      //set this parameter
      currentSFS->paramToWeight.push_back(param);
    }
    else { /* read header */ 
      istringstream headerStream(s);
      vector<int> stateTypeIDs;
      while (! headerStream.eof()) {
	char stateTypeName[STRLEN];
	headerStream.getline(stateTypeName, STRLEN, ',');
	int typeID = getStateTypeID(string(stateTypeName));
	if (typeID == -1) {
	  cerr << "Unrecognized state type name: " << stateTypeName << endl;
	  fatalError("Could not parse DNAKmerFeatures");
	}
        stateTypeIDs.push_back(typeID);
      }
      parameterStream >> k;
      parameterStream >> frameString;
      if (frameString == ".")
	frame = -1;
      else
	frame = atoi(frameString.c_str());
      parameterStream >> speciesName;
      parameterStream >> parameterization;
      sequenceFeatureSets.push_back( new DNAKmerFS(k, frame, parameterization, speciesName) );
      for (int i=0; i<stateTypeIDs.size(); i++) {
	sequenceFeatureSets.back()->stateTypeNames.push_back(stateTypes[stateTypeIDs[i]].name);
	stateTypes[stateTypeIDs[i]].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
      }
      currentSFS = sequenceFeatureSets.back();
    }
  }
}

void SemiCRF::readDNAKmerPairFeatures(ifstream& parameterStream) {
  string s;
  int k;
  string frameString;
  int frame;
  string firstSpeciesName;
  string secondSpeciesName;
  string parameterization;
  SequenceFeatureSet* currentSFS;
  float param;

  //read in the features
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; quit */
      break;
    }
    
    if (sscanf(s.c_str(), "%f", &param)) {
      //set this parameter
      currentSFS->paramToWeight.push_back(param);
    }
    else { /* read header */ 
      istringstream headerStream(s);
      vector<int> stateTypeIDs;
      while (! headerStream.eof()) {
	char stateTypeName[STRLEN];
	headerStream.getline(stateTypeName, STRLEN, ',');
	int typeID = getStateTypeID(string(stateTypeName));
	if (typeID == -1) {
	  cerr << "Unrecognized state type name: " << stateTypeName << endl;
	  fatalError("Could not parse DNAKmerPairFeatures");
	}
        stateTypeIDs.push_back(typeID);
      }
      parameterStream >> k;
      parameterStream >> frameString;
      if (frameString == ".")
	frame = -1;
      else
	frame = atoi(frameString.c_str());
      parameterStream >> firstSpeciesName;
      parameterStream >> secondSpeciesName;
      parameterStream >> parameterization;
      sequenceFeatureSets.push_back( new DNAKmerPairFS(k, frame, parameterization, 
						       firstSpeciesName, secondSpeciesName) );
      for (int i=0; i<stateTypeIDs.size(); i++) {
	sequenceFeatureSets.back()->stateTypeNames.push_back(stateTypes[stateTypeIDs[i]].name);
	stateTypes[stateTypeIDs[i]].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
      }
      currentSFS = sequenceFeatureSets.back();
    }
  }
}

void SemiCRF::readDNAKmerArrayFeatures(ifstream& parameterStream) {
  string s;
  int transitionTypeID;
  int k;
  int length;
  int offset;
  string speciesName;
  string parameterization;
  SequenceFeatureSet* currentSFS;

  //read in the features
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; quit */
      break;
    }
    
    if (getTransitionTypeID(s) != -1) {
      transitionTypeID = getTransitionTypeID(s);
      parameterStream >> k;
      parameterStream >> length;
      parameterStream >> offset;
      parameterStream >> speciesName;
      parameterStream >> parameterization;
      sequenceFeatureSets.push_back(new DNAKmerArrayFS(k, length, offset, parameterization, speciesName));
      transitionTypes[transitionTypeID].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
      currentSFS = transitionTypes[transitionTypeID].sequenceFeatureSets.back();
    }
    else {
      //set this parameter
      currentSFS->paramToWeight.push_back(atof(s.c_str()));
    }
  }
}

void SemiCRF::readDNAKmerPairArrayFeatures(ifstream& parameterStream) {
  string s;
  int transitionTypeID;
  int k;
  int length;
  int offset;
  string firstSpeciesName;
  string secondSpeciesName;
  string parameterization;
  SequenceFeatureSet* currentSFS;

  //read in the features
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; quit */
      break;
    }
    
    if (getTransitionTypeID(s) != -1) {
      transitionTypeID = getTransitionTypeID(s);
      parameterStream >> k;
      parameterStream >> length;
      parameterStream >> offset;
      parameterStream >> firstSpeciesName;
      parameterStream >> secondSpeciesName;
      parameterStream >> parameterization;
      sequenceFeatureSets.push_back(new DNAKmerPairArrayFS(k, length, offset, parameterization, 
							   firstSpeciesName, secondSpeciesName));
      transitionTypes[transitionTypeID].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
      currentSFS = transitionTypes[transitionTypeID].sequenceFeatureSets.back();
    }
    else {
      //set this parameter
      currentSFS->paramToWeight.push_back(atof(s.c_str()));
    }
  }
}

void SemiCRF::readSVMFeatures(ifstream& parameterStream) {
  string s;
  string idString;
  svmfs_t svmfsType;
  kernel_t kernelType;
  int kernelOrder;
  int numBins;
  int length;
  int offset;
  weight_t sampleRate;
  vector<string> seqNames;
  SequenceFeatureSet* currentSFS;
  int id = 0;

  //read in the features
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; quit */
      break;
    }
    else if (getStateTypeID(s) != -1 || getTransitionTypeID(s) != -1) {
      idString = s;
      parameterStream >> s;
      if (s == "Positional")
	svmfsType = POSITIONAL;
      else if (s == "Codon")
	svmfsType = CODON;
      else
	fatalError("Unrecognized SVMFS type in parameter file");
      parameterStream >> s;
      if (s == "Polynomial")
	kernelType = POLYNOMIAL;
      else if (s == "Gaussian")
	kernelType = GAUSSIAN;
      else
	fatalError("Unrecognized kernel type in parameter file");
      parameterStream >> kernelOrder;
      parameterStream >> numBins;
      parameterStream >> length;
      parameterStream >> offset;
      parameterStream >> sampleRate;
      parameterStream >> s;

      //parse sequence names
      seqNames.clear();
      int strPos = 0;
      while (s.find(',', strPos) != string::npos) {
	seqNames.push_back(s.substr(strPos, s.find(',', strPos) - strPos));
	strPos = s.find(',', strPos) + 1;
      }
      seqNames.push_back(s.substr(strPos, s.length() - strPos));

      sequenceFeatureSets.push_back(new SVMFS(id++, svmfsType, kernelType, kernelOrder, numBins, 
					      length, offset, sampleRate, seqNames));
      if (getStateTypeID(idString) != -1) {
	stateTypes[getStateTypeID(idString)].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
	currentSFS = stateTypes[getStateTypeID(idString)].sequenceFeatureSets.back();
      }
      else if (getTransitionTypeID(idString) != -1) {
	transitionTypes[getTransitionTypeID(idString)].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
	currentSFS = transitionTypes[getTransitionTypeID(idString)].sequenceFeatureSets.back();
      }
      else {
	cerr << "Unknown name: " << idString << endl;
	fatalError("Could not find state type or transition type matching SVM feature set name");
      }
    }
    else if (s == "<SVM>") {
      SVMFS* svmfs = (SVMFS*)currentSFS;
      svmfs->readSVM(parameterStream);
    }
    else {
      //set parameter weight and bin limit
      SVMFS* svmfs = (SVMFS*)currentSFS;
      svmfs->paramToWeight.push_back(atof(s.c_str()));
      weight_t binLimit;
      parameterStream >> binLimit;
      svmfs->binLimits.push_back(binLimit);
    }
  }
}

void SemiCRF::readESTPositionFeatures(ifstream& parameterStream) {
  string s;
  int stateTypeID;
  SequenceFeatureSet* currentSFS = NULL;

  //read in the features
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; quit */
      break;
    }
    
    if (getStateTypeID(s) != -1) {
      stateTypeID = getStateTypeID(s);
      sequenceFeatureSets.push_back(new ESTPositionFS());
      stateTypes[stateTypeID].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
      currentSFS = stateTypes[stateTypeID].sequenceFeatureSets.back();
    }
    else if (currentSFS != NULL) {
      //set this parameter
      currentSFS->paramToWeight.push_back(atof(s.c_str()));
    }
  }
}

void SemiCRF::readESTTransitionFeatures(ifstream& parameterStream) {
  string s;
  int transitionTypeID;
  SequenceFeatureSet* currentSFS = NULL;

  //read in the features
  while (parameterStream >> s) {
    if (s[0] == '[') {
      /* we ran into the next header; quit */
      break;
    }
    
    if (getTransitionTypeID(s) != -1) {
      transitionTypeID = getTransitionTypeID(s);
      sequenceFeatureSets.push_back(new ESTTransitionFS());
      transitionTypes[transitionTypeID].sequenceFeatureSets.push_back(sequenceFeatureSets.back());
      currentSFS = transitionTypes[transitionTypeID].sequenceFeatureSets.back();
    }
    else if (currentSFS != NULL) {
      //set this parameter
      currentSFS->paramToWeight.push_back(atof(s.c_str()));
    }
  }
}

void SemiCRF::writeParameters(string filename) {
  ofstream parameterStream(filename.c_str());
  printParameters(parameterStream);
}

void SemiCRF::printParameters(ostream& parameterStream) {
  printStates(parameterStream);
  printTransitions(parameterStream);
  printLengths(parameterStream);
  printMaskingFeatures(parameterStream);
  printDNAKmerFeatures(parameterStream);
  printDNAKmerPairFeatures(parameterStream);
  printDNAKmerArrayFeatures(parameterStream);
  printDNAKmerPairArrayFeatures(parameterStream);
  printSVMFeatures(parameterStream);
  printESTPositionFeatures(parameterStream);
  printESTTransitionFeatures(parameterStream);
}

void SemiCRF::printStates(ostream& parameterStream) {
  parameterStream << "[States]" << endl << endl;
  for (int i=0; i<states.size(); i++) {
    parameterStream << states[i].name << "\t" << stateTypes[states[i].typeID].name <<
      "\t";
    if (states[i].optimize)
      parameterStream << "Yes" << "\t";
    else
      parameterStream << "No" << "\t";
    if (states[i].strand == PLUS)
      parameterStream << "+\t";
    else if (states[i].strand == MINUS)
      parameterStream << "-\t";
    else
      parameterStream << ".\t";
    parameterStream << states[i].startWeight << "\t";
    if (states[i].disallowedKmers.size() == 0)
      parameterStream << ".";
    else
      parameterStream << states[i].disallowedKmers[0];
    for (int j=1; j<states[i].disallowedKmers.size(); j++) {
      parameterStream << "," << states[i].disallowedKmers[j];
    }
    parameterStream << endl;
  }
  parameterStream << endl;
}


void SemiCRF::printTransitions(ostream& parameterStream) {
  parameterStream << "[Transitions]" << endl << endl;
  for (int i=0; i<transitions.size(); i++) {
    parameterStream << states[transitions[i].fromState].name << "\t" <<
      states[transitions[i].toState].name << "\t";

    if (transitions[i].typeID == -1)
      parameterStream << "." << "\t";
    else
      parameterStream << transitionTypes[transitions[i].typeID].name << "\t";

    parameterStream << transitions[i].weight << "\t";

    if(transitions[i].fromEnd.size() == 0)
      parameterStream << ".";
    else
      parameterStream << transitions[i].fromEnd[0];
    for (int j=1; j<transitions[i].fromEnd.size(); j++)
      parameterStream << "," << transitions[i].fromEnd[j];
    parameterStream << "\t";
      
    if(transitions[i].toStart.size() == 0)
      parameterStream << ".";
    else
      parameterStream << transitions[i].toStart[0];
    for (int j=1; j<transitions[i].toStart.size(); j++)
      parameterStream << "," << transitions[i].toStart[j];
    parameterStream << "\t";

    if(transitions[i].disallowedFromEnd.size() == 0)
      parameterStream << ".";
    else
      parameterStream << transitions[i].disallowedFromEnd[0];
    for (int j=1; j<transitions[i].disallowedFromEnd.size(); j++)
      parameterStream << "," << transitions[i].disallowedFromEnd[j];
    parameterStream << "\t";

    if(transitions[i].disallowedToStart.size() == 0)
      parameterStream << ".";
    else
      parameterStream << transitions[i].disallowedToStart[0];
    for (int j=1; j<transitions[i].disallowedToStart.size(); j++)
      parameterStream << "," << transitions[i].disallowedToStart[j];

    parameterStream << endl;
  }
  parameterStream << endl;
}


void SemiCRF::printLengths(ostream& parameterStream) {
  parameterStream << "[Lengths]" << endl << endl;
  for (int i=0; i<stateTypes.size(); i++) {
    if (stateTypes[i].lengthFeatures != NULL) {
      parameterStream << stateTypes[i].name << "\t" << stateTypes[i].lengthFeatures->maxLength <<
	"\t" << stateTypes[i].lengthFeatures->parameterization << endl;
      for (int j=0; j<stateTypes[i].lengthFeatures->paramToWeight.size(); j++)
	parameterStream << stateTypes[i].lengthFeatures->paramToWeight[j] << endl;
      parameterStream << endl;
    }
  }
  parameterStream << endl;
}

void SemiCRF::printMaskingFeatures(ostream& parameterStream) {
  parameterStream << "[MaskingFeatures]" << endl << endl;
  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    if (sequenceFeatureSets[i]->type == MASKING) {
	MaskingFeature* fs = (MaskingFeature*)sequenceFeatureSets[i];
	for (int j=0; j<fs->stateTypeNames.size(); j++)
	  parameterStream << fs->stateTypeNames[j] << (j == fs->stateTypeNames.size()-1 ? "\t" : ",");
	parameterStream << endl;
	for (int k=0; k<fs->paramToWeight.size(); k++)
	  parameterStream << fs->paramToWeight[k] << endl;
	parameterStream << endl;
    }
  }
  parameterStream << endl;
}

void SemiCRF::printDNAKmerFeatures(ostream& parameterStream) {
  parameterStream << "[DNAKmerFeatures]" << endl << endl;
  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    if (sequenceFeatureSets[i]->type == DNAKMER) {
	DNAKmerFS* fs = (DNAKmerFS*)sequenceFeatureSets[i];
	for (int j=0; j<fs->stateTypeNames.size(); j++)
	  parameterStream << fs->stateTypeNames[j] << (j == fs->stateTypeNames.size()-1 ? "\t" : ",");
	parameterStream << fs->k << "\t"
			<< fs->frame << "\t"
			<< fs->seqName << "\t"
			<< fs->parameterization << endl;
	for (int k=0; k<fs->paramToWeight.size(); k++)
	  parameterStream << fs->paramToWeight[k] << endl;
	parameterStream << endl;
    }
  }
  parameterStream << endl;
}

void SemiCRF::printDNAKmerPairFeatures(ostream& parameterStream) {
  parameterStream << "[DNAKmerPairFeatures]" << endl << endl;
  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    if (sequenceFeatureSets[i]->type == DNAKMERPAIR) {
	DNAKmerPairFS* fs = (DNAKmerPairFS*)sequenceFeatureSets[i];
	for (int j=0; j<fs->stateTypeNames.size(); j++)
	  parameterStream << fs->stateTypeNames[j] << (j == fs->stateTypeNames.size()-1 ? "\t" : ",");
	parameterStream << fs->k << "\t"
			<< fs->frame << "\t"
			<< fs->firstSeqName << "\t"
			<< fs->secondSeqName << "\t"
			<< fs->parameterization << endl;
	for (int k=0; k<fs->paramToWeight.size(); k++)
	  parameterStream << fs->paramToWeight[k] << endl;
	parameterStream << endl;
    }
  }
  parameterStream << endl;
}

void SemiCRF::printDNAKmerArrayFeatures(ostream& parameterStream) {
  parameterStream << "[DNAKmerArrayFeatures]" << endl << endl;
  for (int i=0; i<transitionTypes.size(); i++) {
    for (int j=0; j<transitionTypes[i].sequenceFeatureSets.size(); j++) {
      if (transitionTypes[i].sequenceFeatureSets[j]->type == DNAKMERARRAY) {
	DNAKmerArrayFS* fs = (DNAKmerArrayFS*)(transitionTypes[i].sequenceFeatureSets[j]); 
	parameterStream << transitionTypes[i].name << "\t" << 
	  fs->k << "\t" <<
	  fs->length << "\t" <<
	  fs->offset << "\t" <<
	  fs->seqName << "\t" <<
	  fs->parameterization << endl;
	for (int k=0; k<fs->paramToWeight.size(); k++)
	  parameterStream << fs->paramToWeight[k] << endl;
      }
    }
    parameterStream << endl;
  }
  parameterStream << endl;
}

void SemiCRF::printDNAKmerPairArrayFeatures(ostream& parameterStream) {
  parameterStream << "[DNAKmerPairArrayFeatures]" << endl << endl;
  for (int i=0; i<transitionTypes.size(); i++) {
    for (int j=0; j<transitionTypes[i].sequenceFeatureSets.size(); j++) {
      if (transitionTypes[i].sequenceFeatureSets[j]->type == DNAKMERPAIRARRAY) {
	DNAKmerPairArrayFS* fs = (DNAKmerPairArrayFS*)(transitionTypes[i].sequenceFeatureSets[j]); 
	parameterStream << transitionTypes[i].name << "\t" << 
	  fs->k << "\t" <<
	  fs->length << "\t" <<
	  fs->offset << "\t" <<
	  fs->firstSeqName << "\t" <<
	  fs->secondSeqName << "\t" <<
	  fs->parameterization << endl;
	for (int k=0; k<fs->paramToWeight.size(); k++)
	  parameterStream << fs->paramToWeight[k] << endl;
      }
    }
    parameterStream << endl;
  }
  parameterStream << endl;
}

void SemiCRF::printSVMFeatures(ostream& parameterStream) {
  parameterStream << "[SVMFeatures]" << endl << endl;
  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    if (sequenceFeatureSets[i]->type == SVM) {
      //find associated state or transition type
      for (int j=0; j<stateTypes.size(); j++) {
	for (int k=0; k<stateTypes[j].sequenceFeatureSets.size(); k++) {
	  if (stateTypes[j].sequenceFeatureSets[k] == sequenceFeatureSets[i])
	    parameterStream << stateTypes[j].name << "\t";
	}
      }
      for (int j=0; j<transitionTypes.size(); j++) {
	for (int k=0; k<transitionTypes[j].sequenceFeatureSets.size(); k++) {
	  if (transitionTypes[j].sequenceFeatureSets[k] == sequenceFeatureSets[i])
	    parameterStream << transitionTypes[j].name << "\t";
	}
      }
      SVMFS* fs = (SVMFS*)sequenceFeatureSets[i];
      if (fs->svmfsType == POSITIONAL)
	parameterStream << "Positional\t";
      else if (fs->svmfsType == CODON)
	parameterStream << "Codon\t";
      if (fs->kernelType == POLYNOMIAL)
	parameterStream << "Polynomial\t";
      else if (fs->kernelType == GAUSSIAN)
	parameterStream << "Gaussian\t";
      parameterStream << fs->kernelOrder << "\t" <<
	fs->numBins << "\t" <<
	fs->length << "\t" <<
	fs->offset << "\t" <<
	fs->sampleRate << "\t";
      for (int k=0; k<fs->seqNames.size()-1; k++)
	parameterStream << fs->seqNames[k] << ",";
      parameterStream << fs->seqNames.back() << endl;
      for (int k=0; k<fs->paramToWeight.size(); k++)
	parameterStream << fs->paramToWeight[k] << "\t" << fs->binLimits[k] << endl;
      parameterStream << "<SVM>" << endl;
      fs->printSVM(parameterStream);
    }
    parameterStream << endl;
  }
  parameterStream << endl;
}


void SemiCRF::printESTPositionFeatures(ostream& parameterStream) {
  parameterStream << "[ESTPositionFeatures]" << endl << endl;
  for (int i=0; i<stateTypes.size(); i++) {
    for (int j=0; j<stateTypes[i].sequenceFeatureSets.size(); j++) {
      if (stateTypes[i].sequenceFeatureSets[j]->type == ESTPOSITION) {
	ESTPositionFS* fs = (ESTPositionFS*)(stateTypes[i].sequenceFeatureSets[j]); 
	parameterStream << stateTypes[i].name << endl; 
	for (int k=0; k<fs->paramToWeight.size(); k++)
	  parameterStream << fs->paramToWeight[k] << endl;
      }
    }
    parameterStream << endl;
  }
  parameterStream << endl;
}

void SemiCRF::printESTTransitionFeatures(ostream& parameterStream) {
  parameterStream << "[ESTTransitionFeatures]" << endl << endl;
  for (int i=0; i<transitionTypes.size(); i++) {
    for (int j=0; j<transitionTypes[i].sequenceFeatureSets.size(); j++) {
      if (transitionTypes[i].sequenceFeatureSets[j]->type == ESTTRANSITION) {
	ESTTransitionFS* fs = (ESTTransitionFS*)(transitionTypes[i].sequenceFeatureSets[j]); 
	parameterStream << transitionTypes[i].name << endl; 
	for (int k=0; k<fs->paramToWeight.size(); k++)
	  parameterStream << fs->paramToWeight[k] << endl;
      }
    }
    parameterStream << endl;
  }
  parameterStream << endl;
}

void SemiCRF::setAllShortcuts() {
  /* set transition weights */
  for (int i=0; i<transitions.size(); i++) {
    transitions[i].weight = *(parameterMap[transitions[i].globalParamID]);
  }

  /* set valueToWeight arrays for all FeatureSets */
  for (int i=0; i<stateTypes.size(); i++) {
    if (stateTypes[i].lengthFeatures != NULL) {
      stateTypes[i].lengthFeatures->setValueToWeight();
      stateTypes[i].lengthFeatures->setValueToGlobalParamID();
    }
  }
  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    sequenceFeatureSets[i]->setValueToWeight();
    sequenceFeatureSets[i]->setValueToGlobalParamID();
  }
}

//determines which sequences and k values are used for single sequence features
//and pair features
void SemiCRF::getNeededKmersAndSequences(set<int>& dnaKs, set<int>& dnaSeqs, 
					 set<int>& dnaPairKs, set<pair<int,int> >& dnaPairSeqs,
					 set<string>& seqNames) {

  for (int i=0; i<stateTypes.size(); i++) {
    for (int j=0; j<stateTypes[i].sequenceFeatureSets.size(); j++) {
      SequenceFeatureSet* sfs = stateTypes[i].sequenceFeatureSets[j];
      if (sfs->type == DNAKMER) {
	DNAKmerFS* dkfs = (DNAKmerFS*)sfs;
	dnaKs.insert(dkfs->k);
	int id = alignment->getSpeciesID(dkfs->seqName);
	dnaSeqs.insert(id);
	seqNames.insert(dkfs->seqName);
      }
      else if (sfs->type == DNAKMERPAIR) {
	DNAKmerPairFS* dkpfs = (DNAKmerPairFS*)sfs;
	dnaPairKs.insert(dkpfs->k);
	pair<int,int> ids;
	ids.first = alignment->getSpeciesID(dkpfs->firstSeqName);
	ids.second = alignment->getSpeciesID(dkpfs->secondSeqName);
	dnaPairSeqs.insert(ids);
	seqNames.insert(dkpfs->firstSeqName);
	seqNames.insert(dkpfs->secondSeqName);
      }
      else if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	for (int k=0; k<svmfs->seqNames.size(); k++) {
	  seqNames.insert(svmfs->seqNames[k]);
	}
      }
    }
  }
  for (int i=0; i<transitionTypes.size(); i++) {
    for (int j=0; j<transitionTypes[i].sequenceFeatureSets.size(); j++) {
      SequenceFeatureSet* sfs = transitionTypes[i].sequenceFeatureSets[j];
      if (sfs->type == DNAKMERARRAY) {
	DNAKmerArrayFS* dkafs = (DNAKmerArrayFS*)sfs;
	dnaKs.insert(dkafs->k);
	int id = alignment->getSpeciesID(dkafs->seqName);
	dnaSeqs.insert(id);
	seqNames.insert(dkafs->seqName);
      }
      else if (sfs->type == DNAKMERPAIRARRAY) {
	DNAKmerPairArrayFS* dkpafs = (DNAKmerPairArrayFS*)sfs;
	dnaPairKs.insert(dkpafs->k);
	pair<int,int> ids;
	ids.first = alignment->getSpeciesID(dkpafs->firstSeqName);
	ids.second = alignment->getSpeciesID(dkpafs->secondSeqName);
	dnaPairSeqs.insert(ids);
	seqNames.insert(dkpafs->firstSeqName);
	seqNames.insert(dkpafs->secondSeqName);
      }
      else if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	for (int k=0; k<svmfs->seqNames.size(); k++) {
	  seqNames.insert(svmfs->seqNames[k]);
	}
      }
    }
  }
}

void SemiCRF::getNeededSequences(set<string>& seqNames) {

  for (int i=0; i<stateTypes.size(); i++) {
    for (int j=0; j<stateTypes[i].sequenceFeatureSets.size(); j++) {
      SequenceFeatureSet* sfs = stateTypes[i].sequenceFeatureSets[j];
      if (sfs->type == DNAKMER) {
	DNAKmerFS* dkfs = (DNAKmerFS*)sfs;
	seqNames.insert(dkfs->seqName);
      }
      else if (sfs->type == DNAKMERPAIR) {
	DNAKmerPairFS* dkpfs = (DNAKmerPairFS*)sfs;
	seqNames.insert(dkpfs->firstSeqName);
	seqNames.insert(dkpfs->secondSeqName);
      }
      else if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	for (int k=0; k<svmfs->seqNames.size(); k++) {
	  seqNames.insert(svmfs->seqNames[k]);
	}
      }
    }
  }
  for (int i=0; i<transitionTypes.size(); i++) {
    for (int j=0; j<transitionTypes[i].sequenceFeatureSets.size(); j++) {
      SequenceFeatureSet* sfs = transitionTypes[i].sequenceFeatureSets[j];
      if (sfs->type == DNAKMERARRAY) {
	DNAKmerArrayFS* dkafs = (DNAKmerArrayFS*)sfs;
	seqNames.insert(dkafs->seqName);
      }
      else if (sfs->type == DNAKMERPAIRARRAY) {
	DNAKmerPairArrayFS* dkpafs = (DNAKmerPairArrayFS*)sfs;
	seqNames.insert(dkpafs->firstSeqName);
	seqNames.insert(dkpafs->secondSeqName);
      }
      else if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	for (int k=0; k<svmfs->seqNames.size(); k++) {
	  seqNames.insert(svmfs->seqNames[k]);
	}
      }
    }
  }
}


/* precomputes k-mer indices for all the sequences at each position */
void SemiCRF::precomputeKmerIndices() {
  time_t startTime = time(NULL);

  alignment->freeKmerIndices();

  set<int> dnaKs, dnaSeqs, dnaPairKs;
  set<pair<int,int> > dnaPairSeqs;  
  set<string> seqNames;
  getNeededKmersAndSequences(dnaKs, dnaSeqs, dnaPairKs, dnaPairSeqs, seqNames);

  /* first compute the single k-mer indices */
  alignment->maxK = -1;
  for (set<int>::iterator iter = dnaKs.begin(); iter != dnaKs.end(); ++iter) {
    if (*iter > alignment->maxK) alignment->maxK = *iter;
  }
  alignment->maxSeq = -1;
  for (set<int>::iterator iter = dnaSeqs.begin(); iter != dnaSeqs.end(); ++iter) {
    if (*iter > alignment->maxSeq) alignment->maxSeq = *iter;
  }

  alignment->kmerIndices = new int**[alignment->maxK + 1];
  alignment->reverseKmerIndices = new int**[alignment->maxK + 1];
  for (int k=0; k<alignment->maxK+1; k++) {
    alignment->kmerIndices[k] = new int*[2 * (alignment->maxSeq + 1)];
    alignment->reverseKmerIndices[k] = new int*[2 * (alignment->maxSeq + 1)];
    for (int seqID=0; seqID<alignment->maxSeq+1; seqID++) {
      if (dnaKs.find(k) != dnaKs.end() && dnaSeqs.find(seqID) != dnaSeqs.end()) {
	alignment->kmerIndices[k][2 * seqID] = new int[alignment->length];
	alignment->kmerIndices[k][2 * seqID + 1] = new int[alignment->length];
	alignment->reverseKmerIndices[k][2 * seqID] = new int[alignment->length];
	alignment->reverseKmerIndices[k][2 * seqID + 1] = new int[alignment->length];
	for (int pos=0; pos<alignment->length; pos++) {
	  alignment->kmerIndices[k][2 * seqID][pos] = 0;
	  alignment->kmerIndices[k][2 * seqID + 1][pos] = 0;
	}
      }
      else {
	alignment->kmerIndices[k][2 * seqID] = NULL;
	alignment->kmerIndices[k][2 * seqID + 1] = NULL;
      }
    }
  }

  for (int k=0; k<alignment->maxK+1; k++) {
    if (dnaKs.find(k) == dnaKs.end()) continue;
    for (int seqID=0; seqID<alignment->maxSeq+1; seqID++) {
      if (dnaSeqs.find(seqID) == dnaSeqs.end()) continue;
      for (int strand=0; strand<2; strand++) {
	int seqIndex = 2 * seqID + strand;
	const int end = alignment->length - k;
#pragma omp parallel for schedule(static)
	for (int pos=0; pos<end; pos++) {
	  for (int offset=0; offset<k; offset++) {
	    alignment->kmerIndices[k][seqIndex][pos] += DNA_INDEX_COEFF[k - offset - 1] * 
	      DNA_TO_INDEX[alignment->sequenceArray[seqIndex][pos + offset]]; 
	  }
	}
	for (int pos = alignment->length - k; pos<alignment->length; pos++)
	  alignment->kmerIndices[k][seqIndex][pos] = -1;  //out of bounds
	for (int pos=0; pos<alignment->length; pos++) {
	  alignment->reverseKmerIndices[k][seqIndex][pos] = 
	    alignment->kmerIndices[k][seqIndex][alignment->length - 1 - pos];
	}
      }
    }
  }
 
  /* now compute the k-mer pair indices */
  alignment->maxPairK = -1;
  for (set<int>::iterator iter = dnaPairKs.begin(); iter != dnaPairKs.end(); ++iter) {
    if (*iter > alignment->maxPairK) alignment->maxPairK = *iter;
  }
  alignment->maxPairSeq = -1;
  for (set<pair<int,int> >::iterator iter = dnaPairSeqs.begin(); iter != dnaPairSeqs.end(); ++iter) {
    if (iter->first > alignment->maxPairSeq) alignment->maxPairSeq = iter->first;
    if (iter->second > alignment->maxPairSeq) alignment->maxPairSeq = iter->second;
  }

  alignment->kmerPairIndices = new int***[alignment->maxPairK+1];
  alignment->reverseKmerPairIndices = new int***[alignment->maxPairK+1];
  for (int k=0; k<alignment->maxPairK+1; k++) {
    alignment->kmerPairIndices[k] = new int**[2 * (alignment->maxPairSeq + 1)];
    alignment->reverseKmerPairIndices[k] = new int**[2 * (alignment->maxPairSeq + 1)];
    for (int firstSeqID=0; firstSeqID<alignment->maxPairSeq+1; firstSeqID++) {
      alignment->kmerPairIndices[k][2 * firstSeqID] = new int*[2 * (alignment->maxPairSeq + 1)];
      alignment->kmerPairIndices[k][2 * firstSeqID + 1] = new int*[2 * (alignment->maxPairSeq + 1)];
      alignment->reverseKmerPairIndices[k][2 * firstSeqID] = new int*[2 * (alignment->maxPairSeq + 1)];
      alignment->reverseKmerPairIndices[k][2 * firstSeqID + 1] = new int*[2 * (alignment->maxPairSeq + 1)];
      for (int secondSeqID=0; secondSeqID<alignment->maxPairSeq+1; secondSeqID++) {
	pair<int,int> ids(firstSeqID, secondSeqID);
	if (dnaPairKs.find(k) != dnaPairKs.end() && dnaPairSeqs.find(ids) != dnaPairSeqs.end() && firstSeqID < secondSeqID) {
	  alignment->kmerPairIndices[k][2 * firstSeqID][2 * secondSeqID] =  new int[alignment->length];
	  alignment->kmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID + 1] = new int[alignment->length];
	  alignment->reverseKmerPairIndices[k][2 * firstSeqID][2 * secondSeqID] =  new int[alignment->length];
	  alignment->reverseKmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID + 1] = new int[alignment->length];
	  for (int pos=0; pos<alignment->length; pos++) {
	    alignment->kmerPairIndices[k][2 * firstSeqID][2 * secondSeqID][pos] = 0;
	    alignment->kmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID + 1][pos] = 0;
	  }
	}
	else {
	  alignment->kmerPairIndices[k][2 * firstSeqID][2 * secondSeqID] =  NULL;
	  alignment->kmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID + 1] = NULL;
	  alignment->reverseKmerPairIndices[k][2 * firstSeqID][2 * secondSeqID] =  NULL;
	  alignment->reverseKmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID + 1] = NULL;
	}
	alignment->kmerPairIndices[k][2 * firstSeqID][2 * secondSeqID + 1] =  NULL;
	alignment->kmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID] = NULL;
	alignment->reverseKmerPairIndices[k][2 * firstSeqID][2 * secondSeqID + 1] =  NULL;
	alignment->reverseKmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID] = NULL;
      }
    }
  }

  for (int k=0; k<alignment->maxPairK+1; k++) {
    if (dnaPairKs.find(k) == dnaPairKs.end()) continue;
    for (int firstSeqID=0; firstSeqID<alignment->maxPairSeq+1; firstSeqID++) {
      for (int secondSeqID = firstSeqID+1; secondSeqID<alignment->maxPairSeq+1; secondSeqID++) {
	pair<int,int> ids(firstSeqID, secondSeqID);
	if (dnaPairSeqs.find(ids) == dnaPairSeqs.end()) continue;
	for (int strand=0; strand<2; strand++) {
	  int firstSeqIndex = 2 * firstSeqID + strand;
	  int secondSeqIndex = 2 * secondSeqID + strand;
	  const int end = alignment->length - k;
#pragma omp parallel for schedule(static)
	  for (int pos=0; pos<end; pos++) {
	    for (int offset=0; offset<k; offset++) {
	      alignment->kmerPairIndices[k][firstSeqIndex][secondSeqIndex][pos] += 
		DNA_INDEX_COEFF[2 * k - offset - 1] * 
		DNA_TO_INDEX[alignment->sequenceArray[firstSeqIndex][pos + offset]];
	      
	      alignment->kmerPairIndices[k][firstSeqIndex][secondSeqIndex][pos] += 
		DNA_INDEX_COEFF[k - offset - 1] *
		DNA_TO_INDEX[alignment->sequenceArray[secondSeqIndex][pos + offset]];
	    }
	  }
	  for (int pos=alignment->length-k; pos<alignment->length; pos++) {
	    alignment->kmerPairIndices[k][firstSeqIndex][secondSeqIndex][pos] = -1; // out of bounds
	  }
	  for (int pos=0; pos<alignment->length; pos++) {
	    alignment->reverseKmerPairIndices[k][firstSeqIndex][secondSeqIndex][pos] = 
	      alignment->kmerPairIndices[k][firstSeqIndex][secondSeqIndex][alignment->length - 1 - pos];
	  }
	}
      }
    }
  }

  //cerr << "precomputeKmerIndices took " << (time(NULL) - startTime) << " seconds" << endl;
}

/* Precomputes a bunch of stuff that will be needed for dynamic programming.
   This function must be called before beginning DP. */
void SemiCRF::precomputeForDP(Segmentation* fixedPath) {
  time_t startTime = time(NULL);

  /* precompute scores for all state types */
  for (int i=0; i<stateTypes.size(); i++)
    stateTypes[i].scorePositions();

  /* create a matrix that records sums of sequence
     features for each state type with an associated implicit state */
  weight_t** sumsByType = initializeNewWeightMatrix(alignment->length, 2*numberOfImplicitTypes, 0);

  for (int j=0; j<alignment->length; j++) {
    for (int i=0; i<numberOfImplicitTypes; i++) {
      if (stateTypes[i].sequenceFeatureSets.size() == 0) continue;  //nothing to do for this type
      
      sumsByType[j][2*i] += stateTypes[i].plusScores[j];
      sumsByType[j][2*i + 1] += stateTypes[i].minusScores[j];
    }
  }

  
  /* Now set the entries of the stateSFSSums matrix to the appropriate
     entries in the matrix we just computed */
  if (dpMatrix.stateSFSSums == NULL) {
    dpMatrix.stateSFSSums = new weight_t*[alignment->length];
    for (int j=0; j<alignment->length; j++) {
      dpMatrix.stateSFSSums[j] = new weight_t[numberOfImplicitStates];
    }
  }
  
  for (int j=0; j<alignment->length; j++) {
    for (int i=0; i<numberOfImplicitStates; i++) {
      if (fixedPath != NULL && fixedPath->labels[j] != -1 && fixedPath->labels[j] != i)
	dpMatrix.stateSFSSums[j][i] = LOG_ZERO; /* cannot be in this state at this position */
      else {
	if (states[i].strand == PLUS || states[i].strand == NONE) {
	  dpMatrix.stateSFSSums[j][i] = sumsByType[j][2 * states[i].typeID];
	}
	else if (states[i].strand == MINUS) {
	  dpMatrix.stateSFSSums[j][i] = sumsByType[j][2 * states[i].typeID + 1];
	}
      }
    }
  }
    
  /* free the temporary sumsByType matrix */
  for (int j=0; j<alignment->length; j++)
    delete[] sumsByType[j];
  delete[] sumsByType;
  
  setAllowedTransitions();
  
  /* create a matrix that records the sum of SequenceFeatureSets for each transition
     type at each position */
  dpMatrix.transitionSFSSums = new weight_t*[alignment->length];
  for (int j=0; j<alignment->length; j++)
    dpMatrix.transitionSFSSums[j] = new weight_t[2 * transitionTypes.size()];
  
  for (int j=0; j<alignment->length; j++) {
    for (int i=0; i<2 * transitionTypes.size(); i++)
      dpMatrix.transitionSFSSums[j][i] = (weight_t)0;
  }
  
  for (int j=0; j<alignment->length; j++) {
    for (int i=0; i<transitionTypes.size(); i++) {
      if (! dpMatrix.allowedTransitionType[j][i])
	continue;
      
      for (int k=0; k<transitionTypes[i].sequenceFeatureSets.size(); k++) {
	SequenceFeatureSet* sfs = transitionTypes[i].sequenceFeatureSets[k];
	
	if (sfs->type == DNAKMERARRAY) {
	  DNAKmerArrayFS* dkafs = (DNAKmerArrayFS*) sfs;
	  
	  /* first positive strand */
	  pos_t start = j - dkafs->offset;
	  if (start >= 0 && start + dkafs-> length + dkafs->k <= alignment->length) {
	    for (int l=0; l<dkafs->length; l++) {
	      int index = alignment->getDNAKmerIndex(2 * dkafs->seqID, start + l, dkafs->k);
	      index += l * intpow(DNA_CHARS, dkafs->k);
	      dpMatrix.transitionSFSSums[j][2*i] += dkafs->valueToWeight[index];
	    }
	    /* now negative strand */
	    start = alignment->length - j - 1 - dkafs->offset + 1;
	    if (start >= 0 && start + dkafs->length + dkafs->k <= alignment->length) {
	      for (int l=0; l<dkafs->length; l++) {
		int index = alignment->getDNAKmerIndex(2 * dkafs->seqID + 1, start + l, dkafs->k);
		index += l * intpow(DNA_CHARS, dkafs->k);
		dpMatrix.transitionSFSSums[j][2*i + 1] += dkafs->valueToWeight[index];
	      }
	    }
	  }
	}
	
	else if (sfs->type == DNAKMERPAIRARRAY) {
	  DNAKmerPairArrayFS* dkpafs = (DNAKmerPairArrayFS*) sfs;
	  
	  pos_t start = j - dkpafs->offset;
	  if (start >= 0 && start + dkpafs-> length + dkpafs->k <= alignment->length) {
	    for (int l=0; l<dkpafs->length; l++) {
	      int index = alignment->getDNAKmerPairIndex(2 * dkpafs->firstSeqID, 
							 2 * dkpafs->secondSeqID,  
							 start + l, dkpafs->k);
	      index += l * intpow(DNA_CHARS, 2 * dkpafs->k);
	      dpMatrix.transitionSFSSums[j][2*i] += dkpafs->valueToWeight[index];
	    }
	    /* now negative strand */
	    start = alignment->length - j - 1 - dkpafs->offset + 1;
	    if (start >= 0 && start + dkpafs->length + dkpafs->k <= alignment->length) {
	      for (int l=0; l<dkpafs->length; l++) {
		int index = alignment->getDNAKmerPairIndex(2 * dkpafs->firstSeqID + 1, 
							   2 * dkpafs->secondSeqID + 1,
							   start + l, dkpafs->k);
		index += l * intpow(DNA_CHARS, 2 * dkpafs->k);
		dpMatrix.transitionSFSSums[j][2*i + 1] += dkpafs->valueToWeight[index];
	      }
	    }
	  }
	}
	else if (sfs->type == SVM) {
	  SVMFS* svmfs = (SVMFS*) sfs;
	  int plusBin = numericSeq->tags[2 * svmfs->id][j];
	  float plusCoeff = numericSeq->sequenceArray[2 * svmfs->id][j];
	  dpMatrix.transitionSFSSums[j][2*i] += svmfs->valueToWeight[plusBin] * plusCoeff +
	    svmfs->valueToWeight[plusBin+1] * (1 - plusCoeff);

	  int minusBin = numericSeq->tags[2 * svmfs->id + 1][j];
	  float minusCoeff = numericSeq->sequenceArray[2 * svmfs->id + 1][j];
	  dpMatrix.transitionSFSSums[j][2*i + 1] += svmfs->valueToWeight[minusBin] * minusCoeff +
	    svmfs->valueToWeight[minusBin+1] * (1 - minusCoeff);
	}
	else if (sfs->type == ESTTRANSITION) {
	  ESTTransitionFS* esttfs = (ESTTransitionFS*) sfs;
	  int plusIndex = EST_CHARS * EST_TO_INDEX[estSeq->sequence[j-1]] + EST_TO_INDEX[estSeq->sequence[j]]; 
	  int minusIndex = EST_TO_INDEX[estSeq->sequence[j-1]] + EST_CHARS * EST_TO_INDEX[estSeq->sequence[j]];
	  dpMatrix.transitionSFSSums[j][2*i] += esttfs->valueToWeight[plusIndex];
	  dpMatrix.transitionSFSSums[j][2*i + 1] += esttfs->valueToWeight[minusIndex];
	}
      }
    }
  }
  
  for (int i=0; i<stateTypes.size(); i++)
    stateTypes[i].freePositionScores();

  //cerr << "precomputeForDP took " << (time(NULL) - startTime) << " seconds" << endl;
}


void SemiCRF::setAllowedLookup() {
  int maxFromEnd = 0;
  int maxToStart = 0;
  for (int t=0; t<transitions.size(); t++) {
    for (int i=0; i<transitions[t].fromEnd.size(); i++) {
      if (strlen(transitions[t].fromEnd[i]) > maxFromEnd)
	maxFromEnd = strlen(transitions[t].fromEnd[i]);
    }
    for (int i=0; i<transitions[t].disallowedFromEnd.size(); i++) {
      if (strlen(transitions[t].disallowedFromEnd[i]) > maxFromEnd)
	maxFromEnd = strlen(transitions[t].disallowedFromEnd[i]);
    }
    for (int i=0; i<transitions[t].toStart.size(); i++) {
      if (strlen(transitions[t].toStart[i]) > maxToStart)
	maxToStart = strlen(transitions[t].toStart[i]);
    }
    for (int i=0; i<transitions[t].disallowedToStart.size(); i++) {
      if (strlen(transitions[t].disallowedToStart[i]) > maxToStart)
	maxToStart = strlen(transitions[t].disallowedToStart[i]);
    }
  }
  allowedKmerLength = maxFromEnd + maxToStart;
  allowedKmerOffset = maxFromEnd;
  
  allowedLookup = initializeNewBoolMatrix(intpow(DNA_CHARS,allowedKmerLength), 
					  transitions.size(), false);
  
  for (int i=0; i<intpow(DNA_CHARS, allowedKmerLength); i++) {
    char kmer[allowedKmerLength];
    int remainder = i;
    for (int j=0; j<allowedKmerLength; j++) {
      int thisBase = remainder / DNA_INDEX_COEFF[allowedKmerLength - j - 1];
      kmer[j] = INDEX_TO_DNA[thisBase];
      remainder -= thisBase * DNA_INDEX_COEFF[allowedKmerLength - j - 1];
    }
    
    for (int t=0; t<transitions.size(); t++) {
      bool allowedFromEnd = transitions[t].fromEnd.size() == 0 ? true : false;
      bool allowedToStart = transitions[t].toStart.size() == 0 ? true : false;
      for (int j=0; j<transitions[t].fromEnd.size(); j++) {
	if (strncmp(kmer + allowedKmerOffset - strlen(transitions[t].fromEnd[j]), 
		    transitions[t].fromEnd[j], 
		    strlen(transitions[t].fromEnd[j])) == 0)
	    allowedFromEnd = true;
      }
      for (int j=0; j<transitions[t].disallowedFromEnd.size(); j++) {
	if (strncmp(kmer + allowedKmerOffset - strlen(transitions[t].disallowedFromEnd[j]), 
		    transitions[t].disallowedFromEnd[j], 
		    strlen(transitions[t].disallowedFromEnd[j])) == 0)
	    allowedFromEnd = false;
      }
      for (int j=0; j<transitions[t].toStart.size(); j++) {
	if (strncmp(kmer + allowedKmerOffset, transitions[t].toStart[j], 
		    strlen(transitions[t].toStart[j])) == 0)
	    allowedToStart = true;
      }   
      for (int j=0; j<transitions[t].disallowedToStart.size(); j++) {
	if (strncmp(kmer + allowedKmerOffset, transitions[t].disallowedToStart[j], 
		    strlen(transitions[t].disallowedToStart[j])) == 0)
	    allowedToStart = false;
      }
      if (allowedFromEnd && allowedToStart)
	allowedLookup[i][t] = true;
    }

    //disallow any transitions to a state if a disallowedKmer is present
    for (int k=0; k<states.size(); k++) {
      bool allowed = true;

      for (int l=0; l<states[k].disallowedKmers.size(); l++) {
	char* disallowedKmer = states[k].disallowedKmers[l];

	if ( (states[k].strand == PLUS || states[k].strand == NONE) &&
	     strncmp(kmer + allowedKmerOffset + 1 - strlen(disallowedKmer), disallowedKmer, strlen(disallowedKmer)) == 0)
	  allowed = false;
	else if ( states[k].strand == MINUS &&
		  strncmp(kmer + allowedKmerOffset, disallowedKmer, strlen(disallowedKmer)) == 0)
	  allowed = false;
	
	if (! allowed) {
	  for (int t=0; t<states[k].transitionsTo.size(); t++)
	    allowedLookup[i][states[k].transitionsTo[t]] = false;
	}
      }
    }

    /*
    for (int t=0; t<transitions.size(); t++) {
      cerr << i << "\t" << string(kmer,allowedKmerLength) << "\t" << states[transitions[t].fromState].name << "\t"
	   << states[transitions[t].toState].name << "\t" << allowedLookup[i][t] << endl;
    }
    */
  }
}

/* fill in the allowedTransition matrix to determine which transitions are allowed
   at which positions in the sequence */
/* allowedTransition[j][t] = true iff transition t can end at position j */
void SemiCRF::setAllowedTransitions() {
  time_t startTime = time(NULL);

  if (dpMatrix.allowedTransition != NULL)
    fatalError("Tried to set allowedTransition, but it was not previous freed!");

  dpMatrix.allowedTransition = new bool*[alignment->length];

  int allNs = 0;
  for (int i=0; i<allowedKmerLength; i++)
    allNs += DNA_INDEX_COEFF[i] * DNA_TO_INDEX['N'];

  //set allowed transitions for positions in which we cannot see enough surrounding bases
  //to the same as if we saw all Ns
  for (pos_t j=0; j<allowedKmerOffset; j++)
    dpMatrix.allowedTransition[j] = allowedLookup[allNs];
  for (pos_t j=alignment->length - allowedKmerLength + allowedKmerOffset - 1; j<alignment->length; j++)
    dpMatrix.allowedTransition[j] = allowedLookup[allNs];

  for (pos_t j=allowedKmerOffset; j<alignment->length - allowedKmerLength + allowedKmerOffset - 1; j++) {
    int kmer = 0;
    for (int kmerPos=0; kmerPos<allowedKmerLength; kmerPos++)
      kmer += DNA_INDEX_COEFF[allowedKmerLength - kmerPos - 1] *
	DNA_TO_INDEX[alignment->sequenceArray[0][j - allowedKmerOffset + kmerPos]];
    dpMatrix.allowedTransition[j] = allowedLookup[kmer];
  }

  /* create the allowedTransitionType matrix
     allowedTransitionType[j][t] = true iff at least one transition of type t is allowed to
     end at position j */
  dpMatrix.allowedTransitionType = initializeNewBoolMatrix(alignment->length, transitionTypes.size(), false);
  for (pos_t j=0; j<alignment->length; j++) {
    for (int t=0; t<transitions.size(); t++) {
      if (transitions[t].typeID != -1 && dpMatrix.allowedTransition[j][t])
        dpMatrix.allowedTransitionType[j][transitions[t].typeID] = true;
    }
  }

  //cerr << "setAllowedTransitions took " << (time(NULL) - startTime) << " seconds" << endl;
}


void SemiCRF::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  
  dpMatrix.freeEverything();  //all DP matrix entries are now invalid

  if (this->alignment != NULL) {
    this->alignment->freeUncompressedSequences();
    this->alignment->freeKmerIndices();
  }
  if (this->numericSeq != NULL) {
    this->numericSeq->freeUncompressedSequences();
  }
  if (this->estSeq != NULL) {
    this->estSeq->freeUncompressedSequence();
  }

  this->alignment = alignment;
  alignment->uncompressSequences();
  alignment->createReverseComplements();

  this->numericSeq = numericSeq;
  if (numericSeq != NULL && numericSeq->compressedSequenceArray != NULL)
    numericSeq->uncompressSequences();

  this->estSeq = estSeq;
  if (estSeq != NULL)
    estSeq->uncompressSequence();

  dpMatrix.length = alignment->length;

  precomputeKmerIndices();

  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    sequenceFeatureSets[i]->setSequences(alignment, numericSeq, estSeq);
  }
}

void SemiCRF::setParameterMap() {
  parameterMap.clear();

  /* start weights */
  for (int i=0; i<states.size(); i++) {
    parameterMap.push_back(&(states[i].startWeight));
    states[i].startParamID = parameterMap.size() - 1;
  }

  /* transitions */
  transitionWeights.clear();
  
  for (int i=0; i<transitions.size(); i++)
    transitions[i].globalParamID = -1;

  for (int i=0; i<transitions.size(); i++) {
    if (transitions[i].globalParamID == -1) {
      transitionWeights.push_back(transitions[i].weight);
      transitions[i].globalParamID = parameterMap.size() + transitionWeights.size() - 1;

      //find all transitions whose parameters should be tied to this one
      if (states[transitions[i].fromState].name.substr(0,6) == "Intron" &&
	  states[transitions[i].toState].name.substr(0,6) == "Intron") {
	for (int j=0; j<transitions.size(); j++) {
	  if ( transitions[i].strand == transitions[j].strand &&
	       states[transitions[i].fromState].name.substr(0,7) ==
	       states[transitions[j].fromState].name.substr(0,7) &&
	       states[transitions[i].toState].name.substr(0,7) ==
	       states[transitions[j].toState].name.substr(0,7))
	    transitions[j].globalParamID = parameterMap.size() + transitionWeights.size() - 1;
	}
      }
      else if (states[transitions[i].fromState].name.substr(0,6) == "Intron") {
	for (int j=0; j<transitions.size(); j++) {
	  if (transitions[i].toState == transitions[j].toState &&
	      transitions[i].strand == transitions[j].strand &&
	      states[transitions[i].fromState].name.substr(0,7) ==
	      states[transitions[j].fromState].name.substr(0,7))
	    transitions[j].globalParamID = parameterMap.size() + transitionWeights.size() - 1;
	}
      }
      else if (states[transitions[i].toState].name.substr(0,6) == "Intron") {
	for (int j=0; j<transitions.size(); j++) {
	  if (transitions[i].fromState == transitions[j].fromState &&
	      transitions[i].strand == transitions[j].strand &&
	      states[transitions[i].toState].name.substr(0,7) ==
	      states[transitions[j].toState].name.substr(0,7))
	    transitions[j].globalParamID = parameterMap.size() + transitionWeights.size() - 1;
	}
      }
    }
  }

  for (int i=0; i<transitionWeights.size(); i++)
    parameterMap.push_back(&(transitionWeights[i]));


  for (int i=0; i<stateTypes.size(); i++) {
    if (stateTypes[i].lengthFeatures != NULL) { /* lengths */
      stateTypes[i].lengthFeatures->globalParamID.resize(stateTypes[i].lengthFeatures->paramToWeight.size());
      for (int j=0; j<stateTypes[i].lengthFeatures->paramToWeight.size(); j++) {
	parameterMap.push_back(&(stateTypes[i].lengthFeatures->paramToWeight[j]));
	stateTypes[i].lengthFeatures->globalParamID[j] = parameterMap.size() - 1;
      }
    }
  }
  for (int i=0; i<sequenceFeatureSets.size(); i++) { /* sequence features */
    SequenceFeatureSet* sfs = sequenceFeatureSets[i];
    sfs->globalParamID.resize(sfs->paramToWeight.size());
    for (int k=0; k<sfs->paramToWeight.size(); k++) {
      parameterMap.push_back(&(sfs->paramToWeight[k]));
      sfs->globalParamID[k] = parameterMap.size() - 1;
    }
  }
  setAllShortcuts();
}

void SemiCRF::setParameters(const vector<weight_t>& newParameters) {
  if (newParameters.size() != parameterMap.size())
    fatalError("setParameters: length of newParameters vector does not match number of model parameters");

  dpMatrix.freeEverything();  //all DP matrix entries are now invalid

  /* replace parameters */
  for (int i=0; i<parameterMap.size(); i++)
    *(parameterMap[i]) = newParameters[i];

  /* update FeatureSet shortcuts */
  setAllShortcuts();
}

weight_t SemiCRF::viterbi(Segmentation& traceback) {
  dpMatrix.setSegmentEndLinks();

  /* allocate appropriate matrices */
  if (dpMatrix.viterbi != NULL || dpMatrix.maxState != NULL || dpMatrix.maxStateDuration != NULL)
    fatalError("Tried to allocate and initialize Viterbi matrices, but at least one was not previously freed");

  dpMatrix.viterbi = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  dpMatrix.maxState = new stateid_t*[alignment->length];
  for (pos_t j=0; j<alignment->length; j++)
    dpMatrix.maxState[j] = new stateid_t[states.size()];

  dpMatrix.maxStateDuration = new pos_t*[alignment->length];
  for (pos_t j=0; j<alignment->length; j++)
    dpMatrix.maxStateDuration[j] = new pos_t[states.size()];


  for (pos_t j=0; j<alignment->length; j++) {
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      if (j == 0){
	dpMatrix.viterbi[j][i] = dpMatrix.stateSFSSums[j][i] + states[i].startWeight;
	dpMatrix.maxStateDuration[j][i] = 1;
      } 
      else {
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  stateid_t prevState = transition.fromState;
	  weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  weight_t logProb = dpMatrix.viterbi[j-1][prevState]
	    + transitions[states[i].transitionsTo[k]].weight
	    + dpMatrix.stateSFSSums[j][i] + transitionSFSSum;
	  if (logProb > dpMatrix.viterbi[j][i]) {
	    dpMatrix.viterbi[j][i] = logProb;
	    dpMatrix.maxState[j][i] = transition.fromState;
	    dpMatrix.maxStateDuration[j][i] = 1;
	  }
	}
      }
    }
  }

  /* traceback */
  /* find max state at end of matrix */
  stateid_t maxState;
  pos_t j = alignment->length - 1;
  weight_t max = LOG_ZERO;
  for (int i=0; i<states.size(); i++) {
    if (dpMatrix.viterbi[j][i] > max) {
      max = dpMatrix.viterbi[j][i];
      maxState = i;
    }
  }
  weight_t maxPathWeight = dpMatrix.viterbi[j][maxState];

  while (j >= 0) {
    Segment seg;
    seg.state = maxState;
    seg.end = j;
    seg.start = j - dpMatrix.maxStateDuration[j][maxState] + 1;
    if (j != alignment->length -1 && traceback.segments.back().state == maxState)
      traceback.segments.back().start--;
    else
      traceback.segments.push_back(seg);

    stateid_t oldMaxState = maxState;
    maxState = dpMatrix.maxState[j][maxState];
    j -= dpMatrix.maxStateDuration[j][oldMaxState];
  }

  traceback.length = alignment->length;
  reverse(traceback.segments.begin(), traceback.segments.end());

  return maxPathWeight;
}

weight_t SemiCRF::forward() {
  if (dpMatrix.alpha != NULL)
    fatalError("Forward called but matrix not freed");

  dpMatrix.alpha = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=0; j<alignment->length; j++) {
    for (stateid_t i=0; i<numberOfImplicitStates; i++) {
      if (j == 0){
	dpMatrix.alpha[j][i] = dpMatrix.stateSFSSums[j][i] + states[i].startWeight;
      } 
      else {
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  stateid_t prevState = transition.fromState;
	  weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  logPlusEquals(dpMatrix.alpha[j][i], dpMatrix.alpha[j-1][prevState]
			+ transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);
	}
      }
    }
  }

  return computePartitionFromAlpha();
}

weight_t SemiCRF::backward() {
  if (dpMatrix.beta != NULL)
    fatalError("Backward called but matrix not freed");

  dpMatrix.beta = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (stateid_t i=0; i<states.size(); i++)
    dpMatrix.beta[alignment->length-1][i] = (weight_t)0;

  for (pos_t j=alignment->length-2; j>=0; j--) {
    for (stateid_t i=0; i<states.size(); i++) {
      for (int k=0; k<states[i].transitionsFrom.size(); k++){
	if (! dpMatrix.allowedTransition[j+1][states[i].transitionsFrom[k]])
	  continue;

	TransitionFeature& transition = transitions[states[i].transitionsFrom[k]];
	stateid_t nextState = transition.toState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j+1][2 * transition.typeID + transition.strand];

	if (nextState < numberOfImplicitStates){
	  logPlusEquals (dpMatrix.beta[j][i],
			 dpMatrix.beta[j+1][nextState] 
			 + transition.weight
			 + dpMatrix.stateSFSSums[j+1][nextState]
			 + transitionSFSSum);
	} 
      }
    }
  }

  return computePartitionFromBeta();
}

weight_t SemiCRF::computePartitionFromAlpha() {
  weight_t partition = LOG_ZERO;
  for (stateid_t i = 0; i<states.size(); i++) {
    logPlusEquals(partition, dpMatrix.alpha[alignment->length-1][i]);
  }
  return partition;
}

weight_t SemiCRF::computePartitionFromBeta() {
  weight_t partition = LOG_ZERO;
  for (stateid_t i=0; i<states.size(); i++){
    logPlusEquals (partition, dpMatrix.beta[0][i] + dpMatrix.stateSFSSums[0][i] +
		   states[i].startWeight);
  }
  return partition;
}

void SemiCRF::featureExpectations (vector<weight_t> &counts) {
  time_t startTime = time(NULL);

  counts.clear();
  counts.resize (parameterMap.size(), 0);
  
  weight_t partition = computePartitionFromAlpha();
  int i,j,k,l;
  stateid_t prevState;
  weight_t transitionSFSSum, prob, thisProb;
  strand_t strand;
  vector<weight_t> threadCounts(parameterMap.size(), 0);

#pragma omp parallel private (i,j,k,l,prevState,transitionSFSSum,prob,thisProb,strand) firstprivate(threadCounts)
 {
  #pragma omp for schedule(static)
  for (j=0; j<alignment->length; j++) {
    for (i=0; i<numberOfImplicitStates; i++) {
      prob = 0.0;

      if (j == 0) {
	prob = exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		    dpMatrix.beta[j][i] - partition );

	/* add the expectation for the startWeight feature */
	threadCounts[states[i].startParamID] += prob;
      }
      else {
	for (k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  prevState = transition.fromState;
	  transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  thisProb = exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			  dpMatrix.stateSFSSums[j][i] + 
			  transitionSFSSum + 
			  dpMatrix.beta[j][i] - partition );
	  
	  /* add the expectation for this transition feature */
	  threadCounts[transition.globalParamID] += thisProb;
	  
	  /* add the expectations for all sequence features associated with this transition type */
	  if (transition.typeID != -1) {
	    TransitionType& transitionType = transitionTypes[transition.typeID];
	    for (l=0; l<transitionType.sequenceFeatureSets.size(); l++) {
	      strand = transition.strand;
	      transitionType.sequenceFeatureSets[l]->positionExpectationHelper
		(j, transition.strand, thisProb, threadCounts);
	    }
	  }

	  /* add this transition's contribution to the overall probability */
	  prob += thisProb;
	}
      }
      
      /* add expectations for all sequence features associated with this type */
      for (k=0; k<stateTypes[states[i].typeID].sequenceFeatureSets.size(); k++) {
	stateTypes[states[i].typeID].sequenceFeatureSets[k]->positionExpectationHelper
	  (j, states[i].strand, prob, threadCounts);
      }
    }
  }
  #pragma omp critical
  {
    for (i=0; i<threadCounts.size(); i++)
      counts[i] += threadCounts[i];
  }
 }

 // cerr << "featureExpectations took " << (time(NULL) - startTime) << " seconds" << endl;
}

weight_t SemiCRF::computeAnnotationPosteriors(Segmentation& annotation, vector<weight_t>& posteriors, weight_t lambda) {
  /* assume precomputeForDP, forward, and backward have been called */ 

  posteriors.resize(alignment->length);
  weight_t posteriorSum = 0;
  weight_t partition = computePartitionFromAlpha();
  
  for (int j=0; j<alignment->length; j++) {
    stateid_t stateID = annotation.labels[j];
    posteriors[j] = dpMatrix.alpha[j][stateID] + dpMatrix.beta[j][stateID] - partition;
    if (states[stateID].optimize)
      posteriorSum += lambda * exp(posteriors[j]);
    else
      posteriorSum += exp(posteriors[j]);
  }

  return posteriorSum;
}

weight_t SemiCRF::computeESS(Segmentation& annotation, weight_t lambda) {
  /* assume precomputeForDP, forward, and backward have been called */ 

  weight_t ss = 0;
  weight_t partition = computePartitionFromAlpha();
  
  for (int j=0; j<alignment->length; j++) {
    stateid_t stateID = annotation.labels[j];
  
    if (states[stateID].optimize)
      ss += (lambda + 1.0) * exp(dpMatrix.alpha[j][stateID] + dpMatrix.beta[j][stateID] - partition);

    for (int i=0; i<states.size(); i++) {
      if (! states[i].optimize)
	ss += exp(dpMatrix.alpha[j][i] + dpMatrix.beta[j][i] - partition);
    }
  }
 
  //ss -= alignment->length;

  return ss;
}

weight_t SemiCRF::computeESA(Segmentation& annotation, weight_t lambda) {
  /* assume precomputeForDP, forward, and backward have been called */ 

  weight_t esa = 0;
  weight_t partition = computePartitionFromAlpha();
  weight_t threadESA = 0;

#pragma omp parallel firstprivate(threadESA)
 {
#pragma omp for schedule(static)
  for (int j=1; j<alignment->length; j++) {
    for (int i=0; i<states.size(); i++) {
      for (int k=0; k<states[i].transitionsTo.size(); k++) {
	if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	  continue;

	TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	stateid_t prevState = transition.fromState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	weight_t transitionProb = exp(dpMatrix.alpha[j-1][prevState]
			     + dpMatrix.stateSFSSums[j][i]
			     + transition.weight
			     + transitionSFSSum
			     + dpMatrix.beta[j][i]
			     - partition);

	if (states[i].optimize != states[prevState].optimize)
	  threadESA += transitionProb;

	if (states[annotation.labels[j]].optimize != states[annotation.labels[j-1]].optimize &&
	    ! (annotation.labels[j] == i && annotation.labels[j-1] == prevState) )
	  threadESA += (lambda + 1.0) * transitionProb;
      }
    }
  }
#pragma omp critical
  {
    esa += threadESA;
  }
 }

  return esa;
}

weight_t SemiCRF::computeASA(Segmentation& annotation, weight_t lambda, weight_t gamma, weight_t& feMultiplier) {
  // assume precomputeForDP, forward, and backward have been called

  weight_t asa = 0;
  feMultiplier = 0;  //for use in asaGradient
  weight_t numSignals = 0;
  weight_t partition = computePartitionFromAlpha();
  weight_t threadASA = 0;
  weight_t threadFEM = 0;
  weight_t threadNumSignals = 0;

#pragma omp parallel firstprivate(threadASA, threadFEM, threadNumSignals)
 {
#pragma omp for schedule(static)
  for (int j=1; j<alignment->length; j++) {
    if (states[annotation.labels[j]].optimize != states[annotation.labels[j-1]].optimize)
      threadNumSignals++;
    for (int i=0; i<states.size(); i++) {
      for (int k=0; k<states[i].transitionsTo.size(); k++) {
	if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	  continue;

	TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	stateid_t prevState = transition.fromState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	weight_t transitionProb = exp(dpMatrix.alpha[j-1][prevState]
			     + dpMatrix.stateSFSSums[j][i]
			     + transition.weight
			     + transitionSFSSum
			     + dpMatrix.beta[j][i]
			     - partition);

	if (states[i].optimize != states[prevState].optimize) {
	  threadASA += Q(transitionProb - 0.5, gamma);
	  threadFEM += dQdx(transitionProb - 0.5, gamma) * transitionProb;
	}

	if (states[annotation.labels[j]].optimize != states[annotation.labels[j-1]].optimize &&
	    ! (annotation.labels[j] == i && annotation.labels[j-1] == prevState) ) {
	  threadASA += (lambda + 1.0) * Q(transitionProb - 0.5, gamma);
	  threadFEM += (lambda + 1.0) * dQdx(transitionProb - 0.5, gamma) * transitionProb;
	}	
      }
    }
  }
#pragma omp critical
  {
    asa += threadASA;
    feMultiplier += threadFEM;
    numSignals += threadNumSignals;
  }
 }

 asa -= numSignals * (lambda + 1.0);

 return asa;
}


weight_t SemiCRF::computeESAAlphaStar(Segmentation& annotation, weight_t lambda) {
  if (dpMatrix.alpha != NULL || dpMatrix.alphaStar != NULL)
    fatalError("Tried to allocate and initialize alpha and alphaStar matrices, but one was not previously freed");
  
  dpMatrix.alpha = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);
  dpMatrix.alphaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (stateid_t i=0; i<numberOfImplicitStates; i++)
    dpMatrix.alpha[0][i] = dpMatrix.stateSFSSums[0][i] + states[i].startWeight;

  for (pos_t j=1; j<alignment->length; j++) {
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      for (int k=0; k<states[i].transitionsTo.size(); k++) {
	if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	stateid_t prevState = transition.fromState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	
	logPlusEquals(dpMatrix.alpha[j][i], dpMatrix.alpha[j-1][prevState]
		      + transition.weight + dpMatrix.stateSFSSums[j][i] +
		      transitionSFSSum);

	logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alphaStar[j-1][prevState]
		      + transition.weight + dpMatrix.stateSFSSums[j][i] +
		      transitionSFSSum);

	if (states[i].optimize != states[prevState].optimize)
	  logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			+ transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);

	if (states[annotation.labels[j]].optimize != states[annotation.labels[j-1]].optimize &&
	    ! (annotation.labels[j] == i && annotation.labels[j-1] == prevState) )
	  logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			+ log(lambda+1.0) + transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);
      }
    }
  }

  return computeESANumeratorFromAlphaStar();
}

weight_t SemiCRF::computeESABetaStar(Segmentation& annotation, weight_t lambda) {
  if (dpMatrix.beta != NULL || dpMatrix.betaStar != NULL)
    fatalError("Tried to allocate and initialize beta and betaStar matrix, but one was not previously freed");
  
  dpMatrix.beta = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);
  dpMatrix.betaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (stateid_t i=0; i<states.size(); i++)
    dpMatrix.beta[alignment->length-1][i] = (weight_t)0;

  for (pos_t j=alignment->length-2; j>=0; j--){
    for (stateid_t i=0; i<states.size(); i++){
      for (int k=0; k<states[i].transitionsFrom.size(); k++){
	if (! dpMatrix.allowedTransition[j+1][states[i].transitionsFrom[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsFrom[k]];
	stateid_t nextState = transition.toState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j+1][2 * transition.typeID + transition.strand];

	logPlusEquals (dpMatrix.beta[j][i],
		       dpMatrix.beta[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);

	logPlusEquals (dpMatrix.betaStar[j][i],
		       dpMatrix.betaStar[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);

	if (states[nextState].optimize != states[i].optimize)
	    logPlusEquals (dpMatrix.betaStar[j][i],
			   dpMatrix.beta[j+1][nextState] 
			   + transition.weight
			   + dpMatrix.stateSFSSums[j+1][nextState]
			   + transitionSFSSum);

	if (states[annotation.labels[j+1]].optimize != states[annotation.labels[j]].optimize && 
	    !(annotation.labels[j+1] == nextState && annotation.labels[j] == i) )
	  logPlusEquals (dpMatrix.betaStar[j][i],
			 dpMatrix.beta[j+1][nextState] 
			 + log(lambda+1.0)
			 + transition.weight
			 + dpMatrix.stateSFSSums[j+1][nextState]
			 + transitionSFSSum);
      }
    }
  }

  return computeESANumeratorFromBetaStar(annotation, lambda);
}

/*
weight_t SemiCRF::computeESAAlphaStar(Segmentation& annotation, weight_t lambda) {
  if (dpMatrix.alpha != NULL || dpMatrix.alphaStar != NULL)
    fatalError("Tried to allocate and initialize alpha and alphaStar matrices, but one was not previously freed");
  
  dpMatrix.alpha = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);
  dpMatrix.alphaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (stateid_t i=0; i<numberOfImplicitStates; i++)
    dpMatrix.alpha[0][i] = dpMatrix.stateSFSSums[0][i] + states[i].startWeight;

  for (pos_t j=1; j<alignment->length; j++) {
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      for (int k=0; k<states[i].transitionsTo.size(); k++) {
	if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	stateid_t prevState = transition.fromState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	
	logPlusEquals(dpMatrix.alpha[j][i], dpMatrix.alpha[j-1][prevState]
		      + transition.weight + dpMatrix.stateSFSSums[j][i] +
		      transitionSFSSum);

	logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alphaStar[j-1][prevState]
		      + transition.weight + dpMatrix.stateSFSSums[j][i] +
		      transitionSFSSum);
	if (states[i].optimize == states[prevState].optimize)
	  logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			+ transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);
	else if (annotation.labels[j] == i && annotation.labels[j-1] == prevState)
	  logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			+ log(lambda+1.0) + transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);
      }
    }
  }

  return computeESANumeratorFromAlphaStar();
}

weight_t SemiCRF::computeESABetaStar(Segmentation& annotation, weight_t lambda) {
  if (dpMatrix.beta != NULL || dpMatrix.betaStar != NULL)
    fatalError("Tried to allocate and initialize beta and betaStar matrix, but one was not previously freed");
  
  dpMatrix.beta = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);
  dpMatrix.betaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (stateid_t i=0; i<states.size(); i++)
    dpMatrix.beta[alignment->length-1][i] = (weight_t)0;

  for (pos_t j=alignment->length-2; j>=0; j--){
    for (stateid_t i=0; i<states.size(); i++){
      for (int k=0; k<states[i].transitionsFrom.size(); k++){
	if (! dpMatrix.allowedTransition[j+1][states[i].transitionsFrom[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsFrom[k]];
	stateid_t nextState = transition.toState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j+1][2 * transition.typeID + transition.strand];

	logPlusEquals (dpMatrix.beta[j][i],
		       dpMatrix.beta[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);

	logPlusEquals (dpMatrix.betaStar[j][i],
		       dpMatrix.betaStar[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);
	if (states[nextState].optimize == states[i].optimize)
	    logPlusEquals (dpMatrix.betaStar[j][i],
			   dpMatrix.beta[j+1][nextState] 
			   + transition.weight
			   + dpMatrix.stateSFSSums[j+1][nextState]
			   + transitionSFSSum);
	else if (annotation.labels[j+1] == nextState && annotation.labels[j] == i)
	    logPlusEquals (dpMatrix.betaStar[j][i],
			   dpMatrix.beta[j+1][nextState] 
			   + log(lambda+1.0)
			   + transition.weight
			   + dpMatrix.stateSFSSums[j+1][nextState]
			   + transitionSFSSum);
      }
    }
  }

  return computeESANumeratorFromBetaStar(annotation, lambda);
}
*/

weight_t SemiCRF::computeASAAlphaStar(Segmentation& annotation, weight_t lambda, weight_t gamma) {
  if (dpMatrix.alphaStar != NULL)
    fatalError("Tried to allocate and initialize alphaStar matrices, but it was not previously freed");
  
  dpMatrix.alphaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  weight_t partition = computePartitionFromAlpha();

  for (stateid_t i=0; i<numberOfImplicitStates; i++)
    dpMatrix.alpha[0][i] = dpMatrix.stateSFSSums[0][i] + states[i].startWeight;

  for (pos_t j=1; j<alignment->length; j++) {
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      for (int k=0; k<states[i].transitionsTo.size(); k++) {
	if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	stateid_t prevState = transition.fromState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	
	weight_t transitionProb = exp(dpMatrix.alpha[j-1][prevState]
				      + dpMatrix.stateSFSSums[j][i]
				      + transition.weight
				      + transitionSFSSum
				      + dpMatrix.beta[j][i]
				      - partition);

	logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alphaStar[j-1][prevState]
		      + transition.weight + dpMatrix.stateSFSSums[j][i] +
		      transitionSFSSum);

	if (states[i].optimize != states[prevState].optimize)
	  logPlusEquals(dpMatrix.alphaStar[j][i], 
			dpMatrix.alpha[j-1][prevState]
			+ transition.weight 
			+ dpMatrix.stateSFSSums[j][i] 
			+ transitionSFSSum
			+ log(dQdx(transitionProb - 0.5, gamma)));


	if (states[annotation.labels[j]].optimize != states[annotation.labels[j-1]].optimize &&
	    !(annotation.labels[j] == i && annotation.labels[j-1] == prevState))
	  logPlusEquals(dpMatrix.alphaStar[j][i], 
			dpMatrix.alpha[j-1][prevState]
			+ transition.weight 
			+ dpMatrix.stateSFSSums[j][i] 
			+ transitionSFSSum
			+ log(lambda + 1.0)
			+ log(dQdx(transitionProb - 0.5, gamma)));
      }
    }
  }

  return computeASANumeratorFromAlphaStar();
}

weight_t SemiCRF::computeASABetaStar(Segmentation& annotation, weight_t lambda, weight_t gamma) {
  if (dpMatrix.betaStar != NULL)
    fatalError("Tried to allocate betaStar matrix, but it was not previously freed");
  
  dpMatrix.betaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  weight_t partition = computePartitionFromAlpha();

  for (stateid_t i=0; i<states.size(); i++)
    dpMatrix.beta[alignment->length-1][i] = (weight_t)0;

  for (pos_t j=alignment->length-2; j>=0; j--){
    for (stateid_t i=0; i<states.size(); i++){
      for (int k=0; k<states[i].transitionsFrom.size(); k++){
	if (! dpMatrix.allowedTransition[j+1][states[i].transitionsFrom[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsFrom[k]];
	stateid_t nextState = transition.toState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j+1][2 * transition.typeID + transition.strand];

	weight_t transitionProb = exp(dpMatrix.alpha[j][i]
				      + dpMatrix.stateSFSSums[j+1][nextState]
				      + transition.weight
				      + transitionSFSSum
				      + dpMatrix.beta[j+1][nextState]
				      - partition);

	logPlusEquals (dpMatrix.betaStar[j][i],
		       dpMatrix.betaStar[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);

	if (states[nextState].optimize != states[i].optimize)
	  logPlusEquals (dpMatrix.betaStar[j][i],
			 dpMatrix.beta[j+1][nextState] 
			 + transition.weight
			 + dpMatrix.stateSFSSums[j+1][nextState]
			 + transitionSFSSum
			 + log(dQdx(transitionProb - 0.5, gamma)));

	if (states[annotation.labels[j+1]].optimize != states[annotation.labels[j]].optimize &&
	    !(annotation.labels[j+1] == nextState && annotation.labels[j] == i) )
	  logPlusEquals (dpMatrix.betaStar[j][i],
			 dpMatrix.beta[j+1][nextState]
			 + transition.weight
			 + dpMatrix.stateSFSSums[j+1][nextState]
			 + transitionSFSSum
			 + log(lambda + 1.0)
			 + log(dQdx(transitionProb - 0.5, gamma)));
      }
    }
  }

  return computeASANumeratorFromBetaStar(annotation, lambda, gamma);
}

weight_t SemiCRF::computeAlphaStarPPP(Segmentation& annotation, vector<weight_t>& posteriors) {
  if (dpMatrix.alphaStar != NULL)
    fatalError("Tried to allocate and initialize alphaStar matrix, but it was not previously freed");
  
  dpMatrix.alphaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=0; j<alignment->length; j++) {

    /* first the implicit states */
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      if (j == 0) {
	if (annotation.labels[j] == i)
	  dpMatrix.alphaStar[j][i] = dpMatrix.alpha[j][i] - posteriors[j];
	/* otherwise it stays at log(zero) */
      } 
      else {
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  stateid_t prevState = transition.fromState;
	  weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alphaStar[j-1][prevState]
			+ transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);
	  if (annotation.labels[j] == i)
	    logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			  + transition.weight + dpMatrix.stateSFSSums[j][i] +
			  transitionSFSSum - posteriors[j]);
	}
      }
    }
  }

  return computePosteriorSumNumeratorFromAlphaStar();
}

weight_t SemiCRF::computeBetaStarPPP(Segmentation& annotation, vector<weight_t>& posteriors) {
  if (dpMatrix.betaStar != NULL)
    fatalError("Tried to allocate and initialize betaStar matrix, but it was not previously freed");
  
  dpMatrix.betaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=alignment->length-2; j>=0; j--){
    for (stateid_t i=0; i<states.size(); i++){
      for (int k=0; k<states[i].transitionsFrom.size(); k++){
	if (! dpMatrix.allowedTransition[j+1][states[i].transitionsFrom[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsFrom[k]];
	stateid_t nextState = transition.toState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j+1][2 * transition.typeID + transition.strand];

	logPlusEquals (dpMatrix.betaStar[j][i],
		       dpMatrix.betaStar[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);
	if (annotation.labels[j+1] == nextState) {
	  logPlusEquals (dpMatrix.betaStar[j][i],
			 dpMatrix.beta[j+1][nextState] 
			 + transition.weight
			 + dpMatrix.stateSFSSums[j+1][nextState]
			 + transitionSFSSum
			 - posteriors[j+1]);
	}
      }
    }
  }

  return computePosteriorSumNumeratorFromBetaStar(annotation, 1.0);
}

weight_t SemiCRF::computeAlphaStar(Segmentation& annotation, weight_t lambda) {
  if (dpMatrix.alphaStar != NULL)
    fatalError("Tried to allocate and initialize alphaStar matrix, but it was not previously freed");
  
  dpMatrix.alphaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=0; j<alignment->length; j++) {

    /* first the implicit states */
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      if (j == 0) {
	if (annotation.labels[j] == i) {
	  if (states[i].optimize)
	    dpMatrix.alphaStar[j][i] = log(lambda) + dpMatrix.alpha[j][i];
	  else
	    dpMatrix.alphaStar[j][i] = dpMatrix.alpha[j][i];
	}
	/* otherwise it stays at log(zero) */
      } 
      else {
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  stateid_t prevState = transition.fromState;
	  weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alphaStar[j-1][prevState]
			+ transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);
	  if (annotation.labels[j] == i) {
	    if (states[i].optimize)
	      logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			    + log(lambda) + transition.weight + dpMatrix.stateSFSSums[j][i] +
			    transitionSFSSum);
	    else
	      logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			    + transition.weight + dpMatrix.stateSFSSums[j][i] +
			    transitionSFSSum);
	  }
	}
      }
    }
  }

  return computePosteriorSumNumeratorFromAlphaStar();
}

weight_t SemiCRF::computeBetaStar(Segmentation& annotation, weight_t lambda) {
  if (dpMatrix.betaStar != NULL)
    fatalError("Tried to allocate and initialize betaStar matrix, but it was not previously freed");
  
  dpMatrix.betaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=alignment->length-2; j>=0; j--){
    for (stateid_t i=0; i<states.size(); i++){
      for (int k=0; k<states[i].transitionsFrom.size(); k++){
	if (! dpMatrix.allowedTransition[j+1][states[i].transitionsFrom[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsFrom[k]];
	stateid_t nextState = transition.toState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j+1][2 * transition.typeID + transition.strand];

	logPlusEquals (dpMatrix.betaStar[j][i],
		       dpMatrix.betaStar[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);
	if (annotation.labels[j+1] == nextState) {
	  if (states[nextState].optimize)
	    logPlusEquals (dpMatrix.betaStar[j][i],
			   dpMatrix.beta[j+1][nextState] 
			   + log(lambda)
			   + transition.weight
			   + dpMatrix.stateSFSSums[j+1][nextState]
			   + transitionSFSSum);
	  else
	    logPlusEquals (dpMatrix.betaStar[j][i],
			   dpMatrix.beta[j+1][nextState] 
			   + transition.weight
			   + dpMatrix.stateSFSSums[j+1][nextState]
			   + transitionSFSSum);
	}
      }
    }
  }

  return computePosteriorSumNumeratorFromBetaStar(annotation, lambda);
}


weight_t SemiCRF::computeESSAlphaStar(Segmentation& annotation, weight_t lambda) {
  if (dpMatrix.alphaStar != NULL)
    fatalError("Tried to allocate and initialize alphaStar matrix, but it was not previously freed");
  
  dpMatrix.alphaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=0; j<alignment->length; j++) {

    /* first the implicit states */
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      if (j == 0) {
	if (! states[i].optimize)
	  dpMatrix.alphaStar[j][i] = dpMatrix.alpha[j][i];
	else if (states[i].optimize && annotation.labels[j] == i)
	  dpMatrix.alphaStar[j][i] = log(lambda+1.0) + dpMatrix.alpha[j][i];
	/* otherwise it stays at log(zero) */
      } 
      else {
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  stateid_t prevState = transition.fromState;
	  weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  
	  logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alphaStar[j-1][prevState]
			+ transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);
	  if (! states[i].optimize)
	    logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			  + transition.weight + dpMatrix.stateSFSSums[j][i] +
			  transitionSFSSum);

	  if (annotation.labels[j] == i && states[i].optimize)
	    logPlusEquals(dpMatrix.alphaStar[j][i], dpMatrix.alpha[j-1][prevState]
			  + log(lambda+1.0) + transition.weight + dpMatrix.stateSFSSums[j][i] +
			  transitionSFSSum);
	}
      }
    }
  }

  return computeESSNumeratorFromAlphaStar();
}

weight_t SemiCRF::computeESSBetaStar(Segmentation& annotation, weight_t lambda) {
  if (dpMatrix.betaStar != NULL)
    fatalError("Tried to allocate and initialize betaStar matrix, but it was not previously freed");
  
  dpMatrix.betaStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=alignment->length-2; j>=0; j--){
    for (stateid_t i=0; i<states.size(); i++){
      for (int k=0; k<states[i].transitionsFrom.size(); k++){
	if (! dpMatrix.allowedTransition[j+1][states[i].transitionsFrom[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsFrom[k]];
	stateid_t nextState = transition.toState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j+1][2 * transition.typeID + transition.strand];

	logPlusEquals (dpMatrix.betaStar[j][i],
		       dpMatrix.betaStar[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);
	if (! states[nextState].optimize)
	    logPlusEquals (dpMatrix.betaStar[j][i],
			   dpMatrix.beta[j+1][nextState] 
			   + transition.weight
			   + dpMatrix.stateSFSSums[j+1][nextState]
			   + transitionSFSSum);
	if (annotation.labels[j+1] == nextState && states[nextState].optimize)
	    logPlusEquals (dpMatrix.betaStar[j][i],
			   dpMatrix.beta[j+1][nextState] 
			   + log(lambda+1.0)
			   + transition.weight
			   + dpMatrix.stateSFSSums[j+1][nextState]
			   + transitionSFSSum);
      }
    }
  }

  return computeESSNumeratorFromBetaStar(annotation, lambda);
}

/* returns the gradient of the numerator that appears in the expression for a posterior sum, 
  divided by the partition function */
void SemiCRF::posteriorSumNumeratorGradientPPP(vector<weight_t>& g, Segmentation& annotation, vector<weight_t>& posteriors) {
  /* assume that alpha, beta, alphaStar, and betaStar 
     have been computed properly, and that the segment 
     end links are set correctly */
  g.clear();
  g.resize (parameterMap.size(), 0);
  
  weight_t partition = computePartitionFromAlpha();
  int i,j,k,l;
  stateid_t prevState;
  weight_t transitionSFSSum, prob, thisProb;
  strand_t strand;
  vector<weight_t> threadCounts(parameterMap.size(), 0);

#pragma omp parallel private (i,j,k,l,prevState,transitionSFSSum,prob,thisProb,strand) firstprivate(threadCounts)
 {
  #pragma omp for schedule(static)
  for (j=0; j<alignment->length; j++){
    for (i=0; i<numberOfImplicitStates; i++) {
      prob = 0.0;

      if (j == 0) {
	prob = exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		    dpMatrix.betaStar[j][i] - partition );
	if (annotation.labels[j] == i)
	  prob += 
	          exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		       dpMatrix.beta[j][i] - partition - posteriors[j]);

	/* add the expectation for the startWeight feature */
	threadCounts[states[i].startParamID] += prob;
      }
      else {
	for (k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  prevState = transition.fromState;
	  transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  thisProb = exp( dpMatrix.alphaStar[j-1][prevState] + transition.weight +
			  dpMatrix.stateSFSSums[j][i] + 
			  transitionSFSSum + 
			  dpMatrix.beta[j][i] -
			  partition ) +
	             exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			  dpMatrix.stateSFSSums[j][i] + 
			  transitionSFSSum + 
			  dpMatrix.betaStar[j][i] -
			  partition ); 
	  if (annotation.labels[j] == i) {
	    thisProb += exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			     dpMatrix.stateSFSSums[j][i] + 
			     transitionSFSSum + 
			     dpMatrix.beta[j][i] -
			     partition -
			     posteriors[j]);
	  }
	  
	  /* add the expectation for this transition feature */
	  threadCounts[transition.globalParamID] += thisProb;
	  
	  /* add the expectations for all sequence features associated with this transition type */
	  if (transition.typeID != -1) {
	    TransitionType& transitionType = transitionTypes[transition.typeID];
	    for (l=0; l<transitionType.sequenceFeatureSets.size(); l++) {
	      SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[l];
	      sfs->positionExpectationHelper(j, transition.strand, thisProb, g);
	    }
	  }

	  /* add this transition's contribution to the overall log probability */
	  prob += thisProb;
	}
      }
      
      /* add expectations for all sequence features associated with this type */
      for (k=0; k<stateTypes[states[i].typeID].sequenceFeatureSets.size(); k++) {
	stateTypes[states[i].typeID].sequenceFeatureSets[k]->positionExpectationHelper(j, states[i].strand, prob, g);
      }
    }
  }
  #pragma omp critical
  {
    for (i=0; i<threadCounts.size(); i++)
      g[i] += threadCounts[i];
  }
 }
}

/* returns the gradient of the numerator that appears in the expression for a posterior sum, 
  divided by the partition function */
void SemiCRF::posteriorSumNumeratorGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda) {
  /* assume that alpha, beta, alphaStar, and betaStar 
     have been computed properly, and that the segment 
     end links are set correctly */
  g.clear();
  g.resize (parameterMap.size(), 0);
  
  weight_t partition = computePartitionFromAlpha();
  int i,j,k,l;
  stateid_t prevState;
  weight_t transitionSFSSum, prob, thisProb;
  strand_t strand;
  vector<weight_t> threadCounts(parameterMap.size(), 0);

#pragma omp parallel private (i,j,k,l,prevState,transitionSFSSum,prob,thisProb,strand) firstprivate(threadCounts)
 {
  #pragma omp for schedule(static)
  for (j=0; j<alignment->length; j++){
    for (i=0; i<numberOfImplicitStates; i++) {
      prob = 0.0;

      if (j == 0) {
	prob = exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		    dpMatrix.betaStar[j][i] - partition );
	if (annotation.labels[j] == i)
	  prob += 
	          exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		       dpMatrix.beta[j][i] - partition);

	/* add the expectation for the startWeight feature */
	threadCounts[states[i].startParamID] += prob;
      }
      else {
	for (k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  prevState = transition.fromState;
	  transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  thisProb = exp( dpMatrix.alphaStar[j-1][prevState] + transition.weight +
			  dpMatrix.stateSFSSums[j][i] + 
			  transitionSFSSum + 
			  dpMatrix.beta[j][i] -
			  partition ) +
	             exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			  dpMatrix.stateSFSSums[j][i] + 
			  transitionSFSSum + 
			  dpMatrix.betaStar[j][i] -
			  partition ); 
	  if (annotation.labels[j] == i) {
	    if (states[i].optimize)
	      thisProb += exp( dpMatrix.alpha[j-1][prevState] +
			       log(lambda) +
			       transition.weight +
			       dpMatrix.stateSFSSums[j][i] + 
			       transitionSFSSum + 
			       dpMatrix.beta[j][i] -
			       partition);
	    else
	      thisProb += exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			       dpMatrix.stateSFSSums[j][i] + 
			       transitionSFSSum + 
			       dpMatrix.beta[j][i] -
			       partition);
	  }
	  
	  /* add the expectation for this transition feature */
	  threadCounts[transition.globalParamID] += thisProb;
	  
	  /* add the expectations for all sequence features associated with this transition type */
	  if (transition.typeID != -1) {
	    TransitionType& transitionType = transitionTypes[transition.typeID];
	    for (l=0; l<transitionType.sequenceFeatureSets.size(); l++) {
	      SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[l];
	      sfs->positionExpectationHelper(j, transition.strand, thisProb, threadCounts);
	    }
	  }

	  /* add this transition's contribution to the overall log probability */
	  prob += thisProb;
	}
      }
      
      /* add expectations for all sequence features associated with this type */
      for (k=0; k<stateTypes[states[i].typeID].sequenceFeatureSets.size(); k++) {
	stateTypes[states[i].typeID].sequenceFeatureSets[k]->positionExpectationHelper(j, states[i].strand, prob, threadCounts);
      }
    }
  }
  #pragma omp critical
  {
    for (i=0; i<threadCounts.size(); i++)
      g[i] += threadCounts[i];
  }
 }
}

void SemiCRF::essNumeratorGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda) {
  /* assume that alpha, beta, alphaStar, and betaStar 
     have been computed properly, and that the segment 
     end links are set correctly */
  g.clear();
  g.resize (parameterMap.size(), 0);
  
  weight_t partition = computePartitionFromAlpha();
  int i,j,k,l;
  stateid_t prevState;
  weight_t transitionSFSSum, prob, thisProb;
  strand_t strand;
  vector<weight_t> threadCounts(parameterMap.size(), 0);

#pragma omp parallel private (i,j,k,l,prevState,transitionSFSSum,prob,thisProb,strand) firstprivate(threadCounts)
 {
  #pragma omp for schedule(static)
  for (j=0; j<alignment->length; j++){
    for (i=0; i<numberOfImplicitStates; i++) {
      prob = 0.0;

      if (j == 0) {
	prob = exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		    dpMatrix.betaStar[j][i] - partition );
	if (annotation.labels[j] == i)
	  prob += 
	          exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		       dpMatrix.beta[j][i] - partition);

	/* add the expectation for the startWeight feature */
	threadCounts[states[i].startParamID] += prob;
      }
      else {
	for (k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  prevState = transition.fromState;
	  transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  thisProb = exp( dpMatrix.alphaStar[j-1][prevState] + transition.weight +
			  dpMatrix.stateSFSSums[j][i] + 
			  transitionSFSSum + 
			  dpMatrix.beta[j][i] -
			  partition ) +
	             exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			  dpMatrix.stateSFSSums[j][i] + 
			  transitionSFSSum + 
			  dpMatrix.betaStar[j][i] -
			  partition ); 
	  if (! states[i].optimize)
	    thisProb += exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			     dpMatrix.stateSFSSums[j][i] + 
			     transitionSFSSum + 
			     dpMatrix.beta[j][i] -
			     partition);

	  if (annotation.labels[j] == i && states[i].optimize)
	    thisProb += exp( dpMatrix.alpha[j-1][prevState] +
			     log(lambda+1.0) +
			     transition.weight +
			     dpMatrix.stateSFSSums[j][i] + 
			     transitionSFSSum + 
			     dpMatrix.beta[j][i] -
			     partition);
	  
	  /* add the expectation for this transition feature */
	  threadCounts[transition.globalParamID] += thisProb;
	  
	  /* add the expectations for all sequence features associated with this transition type */
	  if (transition.typeID != -1) {
	    TransitionType& transitionType = transitionTypes[transition.typeID];
	    for (l=0; l<transitionType.sequenceFeatureSets.size(); l++) {
	      SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[l];
	      sfs->positionExpectationHelper(j, transition.strand, thisProb, threadCounts);
	    }
	  }

	  /* add this transition's contribution to the overall log probability */
	  prob += thisProb;
	}
      }
      
      /* add expectations for all sequence features associated with this type */
      for (k=0; k<stateTypes[states[i].typeID].sequenceFeatureSets.size(); k++) {
	stateTypes[states[i].typeID].sequenceFeatureSets[k]->positionExpectationHelper(j, states[i].strand, prob, threadCounts);
      }
    }
  }
  #pragma omp critical
  {
    for (i=0; i<threadCounts.size(); i++)
      g[i] += threadCounts[i];
  }
 }
}

void SemiCRF::esaNumeratorGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda) {
  /* assume that alpha, beta, alphaStar, and betaStar 
     have been computed properly, and that the segment 
     end links are set correctly */
  g.clear();
  g.resize (parameterMap.size(), 0);
  
  weight_t partition = computePartitionFromAlpha();
  int i,j,k,l;
  stateid_t prevState;
  weight_t transitionSFSSum, prob, thisProb;
  strand_t strand;
  vector<weight_t> threadCounts(parameterMap.size(), 0);

#pragma omp parallel private (i,j,k,l,prevState,transitionSFSSum,prob,thisProb,strand) firstprivate(threadCounts)
 {
  #pragma omp for schedule(static)
  for (j=1; j<alignment->length; j++){
    for (i=0; i<numberOfImplicitStates; i++) {
      prob = 0.0;
      for (k=0; k<states[i].transitionsTo.size(); k++) {
	if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	prevState = transition.fromState;
	transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	thisProb = exp( dpMatrix.alphaStar[j-1][prevState] + transition.weight +
			dpMatrix.stateSFSSums[j][i] + 
			transitionSFSSum + 
			dpMatrix.beta[j][i] -
			partition ) +
	  exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
	       dpMatrix.stateSFSSums[j][i] + 
	       transitionSFSSum + 
	       dpMatrix.betaStar[j][i] -
	       partition ); 
	if (states[i].optimize == states[prevState].optimize)
	  thisProb += exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			   dpMatrix.stateSFSSums[j][i] + 
			   transitionSFSSum + 
			   dpMatrix.beta[j][i] -
			   partition);
	
	else if (annotation.labels[j] == i && annotation.labels[j-1] == prevState)
	  thisProb += exp( dpMatrix.alpha[j-1][prevState] +
			   log(lambda+1.0) +
			   transition.weight +
			   dpMatrix.stateSFSSums[j][i] + 
			   transitionSFSSum + 
			   dpMatrix.beta[j][i] -
			   partition);
	
	/* add the expectation for this transition feature */
	threadCounts[transition.globalParamID] += thisProb;
	
	/* add the expectations for all sequence features associated with this transition type */
	if (transition.typeID != -1) {
	  TransitionType& transitionType = transitionTypes[transition.typeID];
	  for (l=0; l<transitionType.sequenceFeatureSets.size(); l++) {
	    SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[l];
	    sfs->positionExpectationHelper(j, transition.strand, thisProb, threadCounts);
	  }
	}
	
	/* add this transition's contribution to the overall log probability */
	prob += thisProb;
      }
    
      /* add expectations for all sequence features associated with this type */
      for (k=0; k<stateTypes[states[i].typeID].sequenceFeatureSets.size(); k++) {
	stateTypes[states[i].typeID].sequenceFeatureSets[k]->positionExpectationHelper(j, states[i].strand, prob, threadCounts);
      }
    }
  }
#pragma omp critical
  {
    for (i=0; i<threadCounts.size(); i++)
      g[i] += threadCounts[i];
  }
 }
}



void SemiCRF::esaGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda, weight_t esa) {
  // assume that alpha, beta, alphaStar, and betaStar 
  // have been computed properly, and that the segment 
  // end links are set correctly
  g.clear();
  g.resize (parameterMap.size(), 0);
  
  weight_t** gradientContributions = initializeNewWeightMatrix(states.size(), alignment->length, 0.0);

  weight_t partition = computePartitionFromAlpha();
  vector<weight_t> threadCounts(parameterMap.size(), 0);

#pragma omp parallel firstprivate(threadCounts)
  {
#pragma omp for schedule(static)
    for (int j=1; j<alignment->length; j++){
      for (int i=0; i<numberOfImplicitStates; i++) {
	weight_t prob = 0.0;
	weight_t featureExpectation = 0.0;
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;
	  
	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  int prevState = transition.fromState;
	  weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  weight_t thisProb = exp( dpMatrix.alphaStar[j-1][prevState] + transition.weight +
				   dpMatrix.stateSFSSums[j][i] + 
				   transitionSFSSum + 
				   dpMatrix.beta[j][i] -
				   partition ) +
	    exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
		 dpMatrix.stateSFSSums[j][i] + 
		 transitionSFSSum + 
		 dpMatrix.betaStar[j][i] -
		 partition );
 
	  if (states[i].optimize != states[prevState].optimize)
	    thisProb += exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
			     dpMatrix.stateSFSSums[j][i] + 
			     transitionSFSSum + 
			     dpMatrix.beta[j][i] -
			     partition);
	  
	  if (states[annotation.labels[j]].optimize != states[annotation.labels[j-1]].optimize &&
		   !(annotation.labels[j] == i && annotation.labels[j-1] == prevState) )
	    thisProb += exp( dpMatrix.alpha[j-1][prevState] +
			     log(lambda+1.0) +
			     transition.weight +
			     dpMatrix.stateSFSSums[j][i] + 
			     transitionSFSSum + 
			     dpMatrix.beta[j][i] -
			     partition);
	  
	  weight_t thisFeatureExpectation = exp( dpMatrix.alpha[j-1][prevState] +
						 transition.weight +
						 dpMatrix.stateSFSSums[j][i] + 
						 transitionSFSSum + 
						 dpMatrix.beta[j][i] -
						 partition);
	  
	  // add the expectation for this transition feature
	  threadCounts[transition.globalParamID] += thisProb - esa*thisFeatureExpectation;
	  
	  // add the expectations for all sequence features associated with this transition type
	  if (transition.typeID != -1) {
	    TransitionType& transitionType = transitionTypes[transition.typeID];
	    for (int l=0; l<transitionType.sequenceFeatureSets.size(); l++) {
	      SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[l];
	      sfs->positionExpectationHelper(j, transition.strand, thisProb - esa*thisFeatureExpectation, threadCounts);
	    }
	  }
	  
	  // add this transition's contribution to the overall log probability
	  prob += thisProb;
	  featureExpectation += thisFeatureExpectation;
	}
	
	gradientContributions[i][j] = prob - esa*featureExpectation;
      }
    }

    //calculate gradient of sequence feature set weights based on gradientContributions matrix
    //we do it here instead of in the previous loops for better locality when accessing the gradient vector
#pragma omp for
    for (int i=0; i<numberOfImplicitStates; i++) {
      const strand_t strand = states[i].strand;
      for (int k=0; k<stateTypes[states[i].typeID].sequenceFeatureSets.size(); k++) {
	SequenceFeatureSet* sfs = stateTypes[states[i].typeID].sequenceFeatureSets[k];
	if (sfs->type == DNAKMERPAIR) {  //special optimization for DNA kmer pair features, which are often the most numerous by far
	  DNAKmerPairFS* dkpfs = (DNAKmerPairFS*)sfs;
	  int* index = dkpfs->kmerPairIndices[strand];
	  const int* indexEnd = index + alignment->length;
	  const vector<int>& valueToGlobalParamID = dkpfs->valueToGlobalParamID;
	  weight_t* gradientContribution = gradientContributions[i];
	  for ( ; index != indexEnd; index++) {
	    if (*index != -1)
	      threadCounts[valueToGlobalParamID[*index]] += *gradientContribution;
	    gradientContribution++;
	  }
	}
	else {
	  for (pos_t j=0; j<alignment->length; j++) {
	    sfs->positionExpectationHelper(j, strand, gradientContributions[i][j], threadCounts);
	  }
	}
      }
    }
  
#pragma omp critical
    {
      for (int i=0; i<threadCounts.size(); i++)
	g[i] += threadCounts[i];
    }
  }
 
  freeWeightMatrix(gradientContributions, states.size(), alignment->length);
}


void SemiCRF::asaGradient(vector<weight_t>& g, Segmentation& annotation, weight_t lambda, weight_t gamma, weight_t feMultiplier) {
  // assume that alpha, beta, alphaStar, and betaStar 
  // have been computed properly, and that the segment 
  // end links are set correctly
  g.clear();
  g.resize (parameterMap.size(), 0);
  
  weight_t** gradientContributions = initializeNewWeightMatrix(states.size(), alignment->length, 0.0);

  weight_t partition = computePartitionFromAlpha();
  vector<weight_t> threadCounts(parameterMap.size(), 0);

#pragma omp parallel firstprivate(threadCounts)
  {
#pragma omp for schedule(static)
    for (int j=1; j<alignment->length; j++){
      for (int i=0; i<numberOfImplicitStates; i++) {
	weight_t prob = 0.0;
	weight_t featureExpectation = 0.0;
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;
	  
	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  int prevState = transition.fromState;
	  weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];

	  weight_t transitionProb = exp( dpMatrix.alpha[j-1][prevState] +
					 transition.weight +
					 dpMatrix.stateSFSSums[j][i] + 
					 transitionSFSSum + 
					 dpMatrix.beta[j][i] -
					 partition);

	  weight_t thisProb = exp( dpMatrix.alphaStar[j-1][prevState] + transition.weight +
				   dpMatrix.stateSFSSums[j][i] + 
				   transitionSFSSum + 
				   dpMatrix.beta[j][i] -
				   partition ) +
	    exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
		 dpMatrix.stateSFSSums[j][i] + 
		 transitionSFSSum + 
		 dpMatrix.betaStar[j][i] -
		 partition ); 

	  if (states[i].optimize != states[prevState].optimize)
	    thisProb += exp( dpMatrix.alpha[j-1][prevState] 
			     + log(dQdx(transitionProb - 0.5, gamma))
			     + transition.weight 
			     + dpMatrix.stateSFSSums[j][i] 
			     + transitionSFSSum
			     + dpMatrix.beta[j][i] 
			     - partition);
	  
	  if (states[annotation.labels[j]].optimize != states[annotation.labels[j-1]].optimize &&
	      !(annotation.labels[j] == i && annotation.labels[j-1] == prevState) )
	    thisProb += exp( dpMatrix.alpha[j-1][prevState] 
			     + log(lambda+1.0)
			     + log(dQdx(transitionProb - 0.5, gamma))
			     + transition.weight 
			     + dpMatrix.stateSFSSums[j][i] 
			     + transitionSFSSum 
			     + dpMatrix.beta[j][i] 
			     - partition);
	  
	  weight_t thisFeatureExpectation = transitionProb;
	  
	  // add the expectation for this transition feature
	  threadCounts[transition.globalParamID] += thisProb - feMultiplier * thisFeatureExpectation;
	  
	  // add the expectations for all sequence features associated with this transition type
	  if (transition.typeID != -1) {
	    TransitionType& transitionType = transitionTypes[transition.typeID];
	    for (int l=0; l<transitionType.sequenceFeatureSets.size(); l++) {
	      SequenceFeatureSet* sfs = transitionType.sequenceFeatureSets[l];
	      sfs->positionExpectationHelper(j, transition.strand, thisProb - feMultiplier*thisFeatureExpectation, threadCounts);
	    }
	  }
	  
	  // add this transition's contribution to the overall log probability
	  prob += thisProb;
	  featureExpectation += thisFeatureExpectation;
	}
	
	gradientContributions[i][j] = prob - feMultiplier * featureExpectation;
      }
    }

    //calculate gradient of sequence feature set weights based on gradientContributions matrix
    //we do it here instead of in the previous loops for better locality when accessing the gradient vector
#pragma omp for
    for (int i=0; i<numberOfImplicitStates; i++) {
      const strand_t strand = states[i].strand;
      for (int k=0; k<stateTypes[states[i].typeID].sequenceFeatureSets.size(); k++) {
	SequenceFeatureSet* sfs = stateTypes[states[i].typeID].sequenceFeatureSets[k];
	if (sfs->type == DNAKMERPAIR) {  //special optimization for DNA kmer pair features, which are often the most numerous by far
	  DNAKmerPairFS* dkpfs = (DNAKmerPairFS*)sfs;
	  int* index = dkpfs->kmerPairIndices[strand];
	  const int* indexEnd = index + alignment->length;
	  const vector<int>& valueToGlobalParamID = dkpfs->valueToGlobalParamID;
	  weight_t* gradientContribution = gradientContributions[i];
	  for ( ; index != indexEnd; index++) {
	    if (*index != -1)
	      threadCounts[valueToGlobalParamID[*index]] += *gradientContribution;
	    gradientContribution++;
	  }
	}
	else {
	  for (pos_t j=0; j<alignment->length; j++) {
	    sfs->positionExpectationHelper(j, strand, gradientContributions[i][j], threadCounts);
	  }
	}
      }
    }
  
#pragma omp critical
    {
      for (int i=0; i<threadCounts.size(); i++)
	g[i] += threadCounts[i];
    }
  }
 
  freeWeightMatrix(gradientContributions, states.size(), alignment->length);
}

weight_t SemiCRF::computePosteriorSumNumeratorFromAlphaStar() {
  weight_t numerator = LOG_ZERO;
  for (stateid_t i = 0; i<states.size(); i++) {
    logPlusEquals(numerator, dpMatrix.alphaStar[alignment->length-1][i]);
  }
  return numerator;
}

weight_t SemiCRF::computePosteriorSumNumeratorFromBetaStar(Segmentation& annotation, weight_t lambda) {
  weight_t numerator = LOG_ZERO;
  for (stateid_t i=0; i<states.size(); i++){
    logPlusEquals (numerator, dpMatrix.betaStar[0][i] + dpMatrix.stateSFSSums[0][i] +
		   states[i].startWeight);
    if (annotation.labels[0] == i) {
      if (states[i].optimize)
	logPlusEquals (numerator, dpMatrix.beta[0][i] + log(lambda) + dpMatrix.stateSFSSums[0][i] +
		       states[i].startWeight);
      else
	logPlusEquals (numerator, dpMatrix.beta[0][i] + dpMatrix.stateSFSSums[0][i] +
		       states[i].startWeight);
    }
  }
  return numerator;
}

weight_t SemiCRF::computeESSNumeratorFromAlphaStar() {
  weight_t numerator = LOG_ZERO;
  for (stateid_t i = 0; i<states.size(); i++) {
    logPlusEquals(numerator, dpMatrix.alphaStar[alignment->length-1][i]);
  }
  return numerator;
}

//this is wrong!
weight_t SemiCRF::computeESSNumeratorFromBetaStar(Segmentation& annotation, weight_t lambda) {
  weight_t numerator = LOG_ZERO;
  for (stateid_t i=0; i<states.size(); i++){
    logPlusEquals (numerator, dpMatrix.betaStar[0][i] + dpMatrix.stateSFSSums[0][i] +
		   states[i].startWeight);
    if (! states[i].optimize)
      logPlusEquals (numerator, dpMatrix.beta[0][i] + dpMatrix.stateSFSSums[0][i] +
		     states[i].startWeight);
    if (states[i].optimize && annotation.labels[0] == i)
      logPlusEquals (numerator, dpMatrix.beta[0][i] + log(lambda+1.0) + dpMatrix.stateSFSSums[0][i] +
		     states[i].startWeight);      
  }
  return numerator;
}

weight_t SemiCRF::computeESANumeratorFromAlphaStar() {
  weight_t numerator = LOG_ZERO;
  for (stateid_t i = 0; i<states.size(); i++) {
    logPlusEquals(numerator, dpMatrix.alphaStar[alignment->length-1][i]);
  }
  return numerator;
}

weight_t SemiCRF::computeESANumeratorFromBetaStar(Segmentation& annotation, weight_t lambda) {
  weight_t numerator = LOG_ZERO;

  return numerator;
}

weight_t SemiCRF::computeASANumeratorFromAlphaStar() {
  weight_t numerator = LOG_ZERO;
  for (stateid_t i = 0; i<states.size(); i++) {
    logPlusEquals(numerator, dpMatrix.alphaStar[alignment->length-1][i]);
  }
  return numerator;
}

weight_t SemiCRF::computeASANumeratorFromBetaStar(Segmentation& annotation, weight_t lambda, weight_t gamma) {
  weight_t numerator = LOG_ZERO;

  return numerator;
}

void SemiCRF::computePositionalPosteriors() {
  time_t startTime = time(NULL);

  /* assume forward and backward have already been called */
  weight_t partition = computePartitionFromAlpha();

  if (dpMatrix.posteriors != NULL)
    fatalError("Tried to allocate and initialize posteriors matrix, but it was not previously freed");
  dpMatrix.posteriors = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  int i,j;

  #pragma omp parallel for private(i,j) schedule(static)
  for (j=0; j<alignment->length; j++) {
    for (i=0; i<states.size(); i++) {
      dpMatrix.posteriors[j][i] = exp(dpMatrix.alpha[j][i] + dpMatrix.beta[j][i] - partition); 
    }
  }

  //  cerr << "computePositionalPosteriors took " << (time(NULL) - startTime) << " seconds" << endl;
}

weight_t SemiCRF::meaDecode(Segmentation& traceback, weight_t kappa) {
  // assume computePositionalPosteriors has already been called

  // allocate appropriate matrices
  if (dpMatrix.meaDecode != NULL || dpMatrix.maxState != NULL || dpMatrix.maxStateDuration != NULL)
    fatalError("Tried to allocate and initialize MEA decoding matrices, but at least one was not previously freed");

  dpMatrix.meaDecode = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  dpMatrix.maxState = new stateid_t*[alignment->length];
  for (pos_t j=0; j<alignment->length; j++)
    dpMatrix.maxState[j] = new stateid_t[states.size()];

  dpMatrix.maxStateDuration = new pos_t*[alignment->length];
  for (pos_t j=0; j<alignment->length; j++)
    dpMatrix.maxStateDuration[j] = new pos_t[states.size()];


  for (pos_t j=0; j<alignment->length; j++) {
    // first the implicit states
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      dpMatrix.meaDecode[j][i] = (weight_t)0;
      if (j == 0){
	dpMatrix.meaDecode[j][i] = dpMatrix.posteriors[j][i];
	dpMatrix.maxStateDuration[j][i] = 1;
      } 
      else {
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  stateid_t prevState = transition.fromState;
	  weight_t multiplier = states[i].optimize ? kappa : (weight_t)1.0;
	  weight_t score = dpMatrix.meaDecode[j-1][prevState] + multiplier * dpMatrix.posteriors[j][i];
	  if (score > dpMatrix.meaDecode[j][i]) {
	    dpMatrix.meaDecode[j][i] = score;
	    dpMatrix.maxState[j][i] = transition.fromState;
	    dpMatrix.maxStateDuration[j][i] = 1;
	  }
	}
      }
    }
  }

  // traceback
  traceback.segments.clear();
  traceback.length = alignment->length;

  // find max state at end of matrix
  stateid_t maxState;
  pos_t j = alignment->length - 1;
  weight_t max = LOG_ZERO;
  for (int i=0; i<states.size(); i++) {
    if (dpMatrix.meaDecode[j][i] > max) {
      max = dpMatrix.meaDecode[j][i];
      maxState = i;
    }
  }
  weight_t maxPathWeight = dpMatrix.meaDecode[j][maxState];

  while (j >= 0) {
    Segment seg;
    seg.state = maxState;
    seg.end = j;
    seg.start = j - dpMatrix.maxStateDuration[j][maxState] + 1;
    if (j != alignment->length -1 && traceback.segments.back().state == maxState)
      traceback.segments.back().start--;
    else
      traceback.segments.push_back(seg);

    stateid_t oldMaxState = maxState;
    maxState = dpMatrix.maxState[j][maxState];
    j -= dpMatrix.maxStateDuration[j][oldMaxState];
  }

  reverse(traceback.segments.begin(), traceback.segments.end());

  return maxPathWeight;
}

/* MEA decode without requirement of legal path
weight_t SemiCRF::meaDecode(Segmentation& traceback) {
  cerr << "meaDecode" << endl;

  weight_t posteriorSum = 0;
  traceback.segments.clear();
  traceback.length = alignment->length;

  for (pos_t j=0; j<alignment->length; j++) {
    weight_t maxPosterior = -1;
    stateid_t maxPosteriorLabel = -1;
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      if (dpMatrix.posteriors[j][i] > maxPosterior) {
	maxPosterior = dpMatrix.posteriors[j][i];
	maxPosteriorLabel = i;
      }
    }
    traceback.addSegment(maxPosteriorLabel, j+1, j+1);
    posteriorSum += maxPosterior;
  }

  return posteriorSum;
}
*/


weight_t SemiCRF::mesaDecode(Segmentation& traceback, weight_t kappa) {
  // assume precomputeForDP, forward, and backward have been called

  // allocate appropriate matrices
  if (dpMatrix.meaDecode != NULL || dpMatrix.maxState != NULL || dpMatrix.maxStateDuration != NULL)
    fatalError("Tried to allocate and initialize MEA decoding matrices, but at least one was not previously freed");

  dpMatrix.meaDecode = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  dpMatrix.maxState = new stateid_t*[alignment->length];
  for (pos_t j=0; j<alignment->length; j++)
    dpMatrix.maxState[j] = new stateid_t[states.size()];

  dpMatrix.maxStateDuration = new pos_t*[alignment->length];
  for (pos_t j=0; j<alignment->length; j++)
    dpMatrix.maxStateDuration[j] = new pos_t[states.size()];

  weight_t partition = computePartitionFromAlpha();

  for (pos_t j=0; j<alignment->length; j++) {
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      if (j == 0){
	dpMatrix.meaDecode[j][i] = 0.0;
	dpMatrix.maxStateDuration[j][i] = 1;
      } 
      else {
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  stateid_t prevState = transition.fromState;
	  weight_t score = dpMatrix.meaDecode[j-1][prevState];
	  if (states[i].optimize != states[prevState].optimize) {
	    weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	      dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	    weight_t transitionProb = exp(dpMatrix.alpha[j-1][prevState]
					  + dpMatrix.stateSFSSums[j][i]
					  + transition.weight
					  + transitionSFSSum
					  + dpMatrix.beta[j][i]
					  - partition);
	    score += kappa * transitionProb - (1.0 - transitionProb);
	  }
	  if (score > dpMatrix.meaDecode[j][i]) {
	    dpMatrix.meaDecode[j][i] = score;
	    dpMatrix.maxState[j][i] = transition.fromState;
	    dpMatrix.maxStateDuration[j][i] = 1;
	  }
	}
      }
    }
  }

  // traceback
  traceback.segments.clear();
  traceback.length = alignment->length;

  // find max state at end of matrix
  stateid_t maxState;
  pos_t j = alignment->length - 1;
  weight_t max = LOG_ZERO;
  for (int i=0; i<states.size(); i++) {
    if (dpMatrix.meaDecode[j][i] > max) {
      max = dpMatrix.meaDecode[j][i];
      maxState = i;
    }
  }
  weight_t maxPathWeight = dpMatrix.meaDecode[j][maxState];

  while (j >= 0) {
    //cerr << "TB, j = " << j << endl;
    Segment seg;
    seg.state = maxState;
    seg.end = j;
    seg.start = j - dpMatrix.maxStateDuration[j][maxState] + 1;
    if (j != alignment->length -1 && traceback.segments.back().state == maxState)
      traceback.segments.back().start--;
    else
      traceback.segments.push_back(seg);

    stateid_t oldMaxState = maxState;
    maxState = dpMatrix.maxState[j][maxState];
    j -= dpMatrix.maxStateDuration[j][oldMaxState];
  }

  reverse(traceback.segments.begin(), traceback.segments.end());

  traceback.setLabels();

  if (! traceback.checkLegality(*this))
    fatalError("MESA traceback was not legal");

  return maxPathWeight;
}


NumericSequence* SemiCRF::createNumericSequence(AlignmentSequence* alignment) {
  time_t startTime = time(NULL);

  int numIDs = 0;
  while (getSVMFSByID(numIDs) != NULL)
    numIDs++;

  NumericSequence* numericSeq = new NumericSequence(numIDs, alignment->length);

  setSequences(alignment, numericSeq, NULL);

  setAllowedTransitions();

  for (int i=0; i<transitionTypes.size(); i++) {
    for (int k=0; k<transitionTypes[i].sequenceFeatureSets.size(); k++) {
      SequenceFeatureSet* sfs = transitionTypes[i].sequenceFeatureSets[k];
      if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	if (svmfs->svm != NULL) {
	  //score all appropriate positions in the alignment and record
	  for (pos_t j=0; j<alignment->length; j++) {
	    //defaults
	    numericSeq->sequenceArray[2 * svmfs->id][j] = 0;
	    numericSeq->sequenceArray[2 * svmfs->id + 1][j] = 0;
	    numericSeq->tags[2 * svmfs->id][j] = 0;
	    numericSeq->tags[2 * svmfs->id + 1][j] = 0;
	  }
	  
	  pos_t j;
	  vector<int> x;
	  int bin;
	  weight_t decisionValue;

#pragma omp parallel private(j, x, bin, decisionValue)
	  {
#pragma omp for schedule(static)
	    for (j=0; j<alignment->length; j++) {
	      if (dpMatrix.allowedTransitionType[j][i]) {
		if (svmfs->buildInputVector(x, j, PLUS)) {
		  decisionValue = svmfs->calcDecisionValue(x);
		  bin = lower_bound(svmfs->binLimits.begin(), svmfs->binLimits.end(), decisionValue) 
		    - svmfs->binLimits.begin() - 1;
		  weight_t coeff;
		  if (bin == 0)
		    bin = 1;
		  else if (bin == svmfs->binLimits.size() - 1)
		    bin = svmfs->binLimits.size() - 2;
		  coeff = (svmfs->binLimits[bin+1] - decisionValue) / (svmfs->binLimits[bin+1] - svmfs->binLimits[bin]);
		  numericSeq->sequenceArray[2 * svmfs->id][j] = (double)coeff;
		  numericSeq->tags[2 * svmfs->id][j] = bin;
		}
		if (svmfs->buildInputVector(x, j, MINUS)) {
		  decisionValue = svmfs->calcDecisionValue(x);
		  bin = lower_bound(svmfs->binLimits.begin(), svmfs->binLimits.end(), decisionValue) 
		    - svmfs->binLimits.begin() - 1;
		  weight_t coeff;
		  if (bin == 0)
		    bin = 1;
		  else if (bin == svmfs->binLimits.size() - 1)
		    bin = svmfs->binLimits.size() - 2;
		  coeff = (svmfs->binLimits[bin+1] - decisionValue) / (svmfs->binLimits[bin+1] - svmfs->binLimits[bin]);
		  numericSeq->sequenceArray[2 * svmfs->id + 1][j] = (double)coeff;
		  numericSeq->tags[2 * svmfs->id + 1][j] = bin;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  for (int i=0; i<stateTypes.size(); i++) {
    for (int j=0; j<stateTypes[i].sequenceFeatureSets.size(); j++) {
      SequenceFeatureSet* sfs = stateTypes[i].sequenceFeatureSets[j];
      if (sfs->type == SVM) {
	SVMFS* svmfs = (SVMFS*)sfs;
	if (svmfs->svm != NULL) {
	  //score ALL positions in the alignment and record
	  for (pos_t j=0; j<alignment->length; j++) {
	    //default score
	    numericSeq->sequenceArray[2 * svmfs->id][j] = 0;
	    numericSeq->sequenceArray[2 * svmfs->id + 1][j] = 0;
	  }
	  
	  vector<int> x;
	  for (pos_t j=0; j<alignment->length; j++) {
	    if (svmfs->buildInputVector(x, j, PLUS)) {
	      weight_t decisionValue = svmfs->calcDecisionValue(x);
	      int bin = lower_bound(svmfs->binLimits.begin(), svmfs->binLimits.end(), decisionValue) - svmfs->binLimits.begin() - 1;
	      numericSeq->sequenceArray[2 * svmfs->id][j] = (double)decisionValue;
	      numericSeq->tags[2 * svmfs->id][j] = bin;
	    }
	    if (svmfs->buildInputVector(x, j, MINUS)) {
	      weight_t decisionValue = svmfs->calcDecisionValue(x);
	      int bin = lower_bound(svmfs->binLimits.begin(), svmfs->binLimits.end(), decisionValue) - svmfs->binLimits.begin() - 1;
	      numericSeq->sequenceArray[2 * svmfs->id + 1][j] = (double)decisionValue;	      
	      numericSeq->tags[2 * svmfs->id + 1][j] = bin;
	    }
	  }
	}
      }
    }
  }

  numericSeq->compressSequences();
  numericSeq->freeUncompressedSequences();

  //cerr << "createNumericSequence took " << (time(NULL) - startTime) << " seconds" << endl;

  return numericSeq;
}



void SemiCRF::computeAlphaStarStar(Segmentation& annotation, Segmentation& correctAnnotation, 
				   vector<weight_t>& posteriorDifference, weight_t gamma, weight_t lambda) {
  time_t startTime = time(NULL);

  if (dpMatrix.alphaStarStar != NULL)
    fatalError("Tried to allocate and initialize alphaStarStar matrix, but it was not previously freed");
  
  dpMatrix.alphaStarStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=0; j<alignment->length; j++) {

    /* first the implicit states */
    for (stateid_t i=0; i<numberOfImplicitStates; i++){
      if (j == 0) {
	if (annotation.labels[j] == i)
	  dpMatrix.alphaStarStar[j][i] = dpMatrix.alpha[j][i] + log(dQdx(posteriorDifference[j], gamma));
	/* otherwise it stays at log(zero) */
      } 
      else {
	for (int k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  stateid_t prevState = transition.fromState;
	  weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  logPlusEquals(dpMatrix.alphaStarStar[j][i], dpMatrix.alphaStarStar[j-1][prevState]
			+ transition.weight + dpMatrix.stateSFSSums[j][i] +
			transitionSFSSum);
	  if (annotation.labels[j] == i) {
	    weight_t multiplier = states[correctAnnotation.labels[j]].optimize ? lambda : (weight_t)1.0;
	    logPlusEquals(dpMatrix.alphaStarStar[j][i], dpMatrix.alpha[j-1][prevState]
			  + log(multiplier) + transition.weight + dpMatrix.stateSFSSums[j][i] +
			  transitionSFSSum + log(dQdx(posteriorDifference[j], gamma)));
	  }
	}
      }
    }
  }

  //  cerr << "computeAlphaStarStar took " << (time(NULL) - startTime) << " seconds" << endl;
}

void SemiCRF::computeBetaStarStar(Segmentation& annotation, Segmentation& correctAnnotation, 
				  vector<weight_t>& posteriorDifference, weight_t gamma, weight_t lambda) {
  time_t startTime = time(NULL);

  if (dpMatrix.betaStarStar != NULL)
    fatalError("Tried to allocate and initialize betaStarStar matrix, but it was not previously freed");
  
  dpMatrix.betaStarStar = initializeNewWeightMatrix(alignment->length, states.size(), LOG_ZERO);

  for (pos_t j=alignment->length-2; j>=0; j--){
    for (stateid_t i=0; i<states.size(); i++){
      for (int k=0; k<states[i].transitionsFrom.size(); k++){
	if (! dpMatrix.allowedTransition[j+1][states[i].transitionsFrom[k]])
	  continue;
	
	TransitionFeature& transition = transitions[states[i].transitionsFrom[k]];
	stateid_t nextState = transition.toState;
	weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	  dpMatrix.transitionSFSSums[j+1][2 * transition.typeID + transition.strand];

	logPlusEquals (dpMatrix.betaStarStar[j][i],
		       dpMatrix.betaStarStar[j+1][nextState] 
		       + transition.weight
		       + dpMatrix.stateSFSSums[j+1][nextState]
		       + transitionSFSSum);
	if (annotation.labels[j+1] == nextState) {
	  weight_t multiplier = states[correctAnnotation.labels[j+1]].optimize ? lambda : (weight_t)1.0;
	  logPlusEquals (dpMatrix.betaStarStar[j][i],
			 dpMatrix.beta[j+1][nextState]
			 + log(multiplier)
			 + transition.weight
			 + dpMatrix.stateSFSSums[j+1][nextState]
			 + transitionSFSSum
			 + log(dQdx(posteriorDifference[j+1], gamma)));
	}
      }
    }
  }

  //  cerr << "computeBetaStarStar took " << (time(NULL) - startTime) << " seconds" << endl;
}

void SemiCRF::computePosteriorDifference(Segmentation& annotation, vector<weight_t>& posteriorDifference) {
  posteriorDifference.resize(alignment->length);

  int i,j;
  weight_t correctPosterior, maxIncorrectPosterior;

  for (j=0; j<alignment->length; j++) {
    correctPosterior = dpMatrix.posteriors[j][annotation.labels[j]];
    maxIncorrectPosterior = -1;
    for (i=0; i<states.size(); i++) {
      if (i != annotation.labels[j] && dpMatrix.posteriors[j][i] > maxIncorrectPosterior)
	maxIncorrectPosterior = dpMatrix.posteriors[j][i];
    }  
    posteriorDifference[j] = correctPosterior - maxIncorrectPosterior;
  }
}

void SemiCRF::getMaxPosteriorIncorrectLabels(Segmentation& annotation, Segmentation& MPIL) {
  MPIL.segments.clear();
  MPIL.length = alignment->length;

  for (pos_t j=0; j<alignment->length; j++) {
    weight_t correctPosterior = dpMatrix.posteriors[j][annotation.labels[j]];
    stateid_t maxIncorrectPosteriorLabel = -1;
    weight_t maxIncorrectPosterior = -1;
    for (stateid_t i=0; i<states.size(); i++) {
      if (i != annotation.labels[j] && dpMatrix.posteriors[j][i] > maxIncorrectPosterior) {
	maxIncorrectPosteriorLabel = i;
	maxIncorrectPosterior = dpMatrix.posteriors[j][i];
      }
    }  
    MPIL.addSegment(maxIncorrectPosteriorLabel, j+1, j+1);  //addSegment wants 1-based coordinates
  }
  MPIL.setLabels();
}

weight_t SemiCRF::AEasyGradTerm(Segmentation& annotation, Segmentation& correctAnnotation,
				vector<weight_t>& posteriorDifference, weight_t gamma, weight_t lambda) {
  weight_t term = 0;
  for (pos_t j=0; j<alignment->length; j++) {
    if (states[correctAnnotation.labels[j]].optimize)
      term += dQdx(posteriorDifference[j], gamma) * lambda * dpMatrix.posteriors[j][annotation.labels[j]];
    else
      term += dQdx(posteriorDifference[j], gamma) * dpMatrix.posteriors[j][annotation.labels[j]];
  } 

  return term;
}


/* returns the gradient of the numerator that appears in the expression for a posterior sum, 
  divided by the partition function */
void SemiCRF::AHardGradTerm(Segmentation& annotation, Segmentation& correctAnnotation, 
			    vector<weight_t>& posteriorDifference, weight_t gamma, 
			    weight_t lambda, vector<weight_t>& g) {
  time_t startTime = time(NULL);

  g.clear();
  g.resize (parameterMap.size(), 0);
  
  weight_t partition = computePartitionFromAlpha();

  vector<weight_t> threadG(parameterMap.size(), 0);

  int i,j,k,l;
  stateid_t prevState;
  weight_t transitionSFSSum, prob, thisProb;
  strand_t strand;

#pragma omp parallel private(i,j,k,l,prevState,transitionSFSSum,prob,thisProb,strand) firstprivate(threadG)
 {
  #pragma omp for schedule(static)
  for (j=0; j<alignment->length; j++){
    for (i=0; i<numberOfImplicitStates; i++) {
      prob = 0.0;

      if (j == 0) {
	prob = exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		    dpMatrix.betaStarStar[j][i] - partition );
	if (annotation.labels[j] == i)
	  prob += exp( states[i].startWeight + dpMatrix.stateSFSSums[j][i] + 
		    dpMatrix.beta[j][i] - partition );

	/* add the expectation for the startWeight feature */
	threadG[states[i].startParamID] += prob;
      }
      else {
	for (k=0; k<states[i].transitionsTo.size(); k++) {
	  if (! dpMatrix.allowedTransition[j][states[i].transitionsTo[k]])
	    continue;

	  TransitionFeature& transition = transitions[states[i].transitionsTo[k]];
	  prevState = transition.fromState;
	  transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
	    dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
	  thisProb = exp( dpMatrix.alphaStarStar[j-1][prevState] + transition.weight +
				   dpMatrix.stateSFSSums[j][i] + 
				   transitionSFSSum + 
				   dpMatrix.beta[j][i] -
				   partition ) +
	             exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
				   dpMatrix.stateSFSSums[j][i] + 
				   transitionSFSSum + 
				   dpMatrix.betaStarStar[j][i] -
				   partition ); 
	  if (annotation.labels[j] == i) {
	    weight_t multiplier = states[correctAnnotation.labels[j]].optimize ? lambda : (weight_t)1.0;
	    thisProb += 
	      dQdx(posteriorDifference[j], gamma) *
	      multiplier *
	      exp( dpMatrix.alpha[j-1][prevState] + transition.weight +
		   dpMatrix.stateSFSSums[j][i] + 
		   transitionSFSSum + 
		   dpMatrix.beta[j][i] -
		   partition );
	  }
	  
	  /* add the expectation for this transition feature */
	  threadG[transition.globalParamID] += thisProb;
	  
	  /* add the expectations for all sequence features associated with this transition type */
	  if (transition.typeID != -1) {
	    TransitionType& transitionType = transitionTypes[transition.typeID];
	    for (l=0; l<transitionType.sequenceFeatureSets.size(); l++) {
	      transitionType.sequenceFeatureSets[l]->
		positionExpectationHelper(j, transition.strand, thisProb, threadG);
	    }
	  }

	  /* add this transition's contribution to the overall log probability */
	  prob += thisProb;
	}
      }
      
      /* add expectations for all sequence features associated with this type */
      for (k=0; k<stateTypes[states[i].typeID].sequenceFeatureSets.size(); k++) {
	stateTypes[states[i].typeID].sequenceFeatureSets[k]->positionExpectationHelper
	  (j, states[i].strand, prob, threadG);
      }
    }
  }
  #pragma omp critical
  {
    for (i=0; i<threadG.size(); i++)
      g[i] += threadG[i];
  }
 }

 //  cerr << "AHardGradTerm took " << (time(NULL) - startTime) << " seconds" << endl;
}
