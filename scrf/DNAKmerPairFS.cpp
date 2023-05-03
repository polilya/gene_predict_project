#include <math.h>
#include "DNAKmerPairFS.h"
#include "KmerIterator.h"


DNAKmerPairFS::DNAKmerPairFS(int k, int frame, string parameterization, 
			     string firstSeqName, string secondSeqName) {
  this->type = DNAKMERPAIR;
  this->k = k; 
  this->frame = frame;
  this->parameterization = parameterization;
  this->firstSeqName = firstSeqName;
  this->secondSeqName = secondSeqName;
  firstSeqID = -1;  //unknown as yet
  secondSeqID = -1;

  int numValues = intpow(DNA_CHARS, 2*k);
  valueToParam.resize(numValues);
  valueToWeight.resize(numValues);
  valueToGlobalParamID.resize(numValues);

  if (parameterization == "Standard") {
    numParams = intpow(4,2*k) + (intpow(2,k)-1) + 3 + 1;

    KmerIterator kmerIter(2*k, "DNA");
    for (kmerIter.init(); ! kmerIter.finished; kmerIter.next()) {
      int firstKmerParam = kmerToParameter(kmerIter.kmer, k);
      int secondKmerParam = kmerToParameter(kmerIter.kmer + k, k);

      if (firstKmerParam >= intpow(4,k))
	valueToParam[kmerIter.index] = numParams - 1;
      else if (secondKmerParam < intpow(4,k))
	valueToParam[kmerIter.index] = intpow(4,k) * firstKmerParam + secondKmerParam;
      else
	valueToParam[kmerIter.index] = intpow(4,2*k) + (secondKmerParam - intpow(4,k));
      //cerr << kmerIter.kmer << "\t" << kmerIter.index << "\t" << firstKmerParam << "\t" 
      //	   << secondKmerParam << "\t" << valueToParam[kmerIter.index] << endl;
    }
  }
  else if (parameterization == "General") {
    numParams = numValues;
    for (int i=0; i<numValues; i++)
      valueToParam[i] = i;
  }
  else
    fatalError("Unknown parameterization for DNAKmerPairFS");
}

//gets the parameter number of an unpaired k-mer
int DNAKmerPairFS::kmerToParameter(char* kmer, int k) {
  int numKmerParams = intpow(4,k) + (intpow(2,k)-1) + 3;
  int parameter = 0;

  int* indices = new int[k];
  for (int i=0; i<k; i++)
    indices[i] = DNA_TO_INDEX[kmer[i]];

  bool containsGaps = false;
  bool containsUnaligned = false;
  bool containsNs = false;
  for (int i=0; i<k; i++) {
    if (kmer[i] == '_')
      containsGaps = true;
    if (kmer[i] == '.')
      containsUnaligned = true;
    if (kmer[i] == 'N')
      containsNs = true;
  }
  if (!containsGaps && !containsUnaligned && !containsNs) {
    //parameter based on DNA k-mer
    parameter = 0;
    for (int i=0; i<k; i++)
      parameter += intpow(4,i) * indices[k - 1 - i];
  }
  else if (!containsUnaligned && !containsNs) {
    //parameter based on gap pattern
    parameter = intpow(4,k) - 1;
    for (int i=0; i<k; i++) {
      if (kmer[k - 1 - i] == '_')
	parameter += intpow(2,i);
    }
  }
  else if (!containsNs) {
    bool allUnaligned = true;
    for (int i=0; i<k; i++) {
      if (kmer[i] != '.')
	allUnaligned = false;
    }
    if (allUnaligned)
      parameter = numKmerParams - 1;
    else {
      bool alignmentStart = false;
      bool alignmentEnd = false;

      for (int i=1; i<k && !alignmentStart; i++) {  //check for start of alignment at each possible position
	alignmentStart = true;

	if (kmer[i] != 'A' && kmer[i] != 'C' &&kmer[i] != 'G' && kmer[i] != 'T')
	  alignmentStart = false;

	for (int j=i-1; j>=0; j--) {
	  if (kmer[j] != '.')
	    alignmentStart = false;
	}
	for (int j=i+1; j<k; j++) {
	  if (kmer[j] == '.')
	    alignmentStart = false;
	}
      }
      for (int i=0; i<k-1 && !alignmentStart && !alignmentEnd; i++) {
        alignmentEnd = true;

	if (kmer[i] != 'A' && kmer[i] != 'C' && kmer[i] != 'G' && kmer[i] != 'T')
	  alignmentEnd = false;

	for (int j=i+1; j<k; j++) {
	  if (kmer[j] != '.')
	    alignmentEnd = false;
	} 

	for (int j=i-1; j>=0; j--) {
	  if (kmer[j] == '.')
	    alignmentEnd = false;
	}
      }
      if (alignmentStart)
	parameter = numKmerParams - 3;
      else if (alignmentEnd)
	parameter = numKmerParams - 2;
      else
	parameter = numKmerParams;
    }
  }
  else
    parameter = numKmerParams;

  delete[] indices;      
  return parameter;
}

void DNAKmerPairFS::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  firstSeqID = alignment->getSpeciesID(firstSeqName);
  secondSeqID = alignment->getSpeciesID(secondSeqName);
  kmerPairIndices[PLUS] = kmerPairIndices[NONE] = alignment->kmerPairIndices[k][2 * firstSeqID][2 * secondSeqID];
  kmerPairIndices[MINUS] = alignment->reverseKmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID + 1];
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;
}

void DNAKmerPairFS::runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				     vector< vector<weight_t> >& runningSumMinus) {
  int plusFirstRow = 2 * firstSeqID;
  int minusFirstRow = 2 * firstSeqID + 1;
  int plusSecondRow = 2 * secondSeqID;
  int minusSecondRow = 2 * secondSeqID + 1;

  vector<weight_t> ourRunningSumPlus(3, 0);
  vector<weight_t> ourRunningSumMinus(3, 0);

  for (pos_t j=0; j<alignment->length - k + 1; j++) {
    weight_t plusContribution = valueToWeight[kmerPairIndices[PLUS][j]];
    weight_t minusContribution = valueToWeight[kmerPairIndices[MINUS][j]];

    if (frame == -1) {
      runningSumPlus[0][j] += plusContribution + ourRunningSumPlus[0];
      runningSumMinus[0][j] += minusContribution + ourRunningSumMinus[0];
      ourRunningSumPlus[0] += plusContribution;
      ourRunningSumMinus[0] += minusContribution;
    }
    else {
      for (int i=0; i<3; i++) {
	if ((i + (j % 3)) % 3 == frame) {
	  runningSumPlus[i][j] += plusContribution + ourRunningSumPlus[i];
	  runningSumMinus[i][j] += minusContribution + ourRunningSumMinus[i];
	  ourRunningSumPlus[i] += plusContribution;
	  ourRunningSumMinus[i] += minusContribution;
	}
	else { /* out of frame; contribution is zero */
	  runningSumPlus[i][j] += ourRunningSumPlus[i];
	  runningSumMinus[i][j] += ourRunningSumMinus[i];
	}
      }
    }
  }

  /* contribution is zero for out of bounds positions */
  for (pos_t j=alignment->length - k + 1; j<alignment->length; j++) {
    for (int i=0; i<runningSumPlus.size(); i++) {
      runningSumPlus[i][j] += ourRunningSumPlus[i];
      runningSumMinus[i][j] += ourRunningSumMinus[i];
    }
  }
}


void DNAKmerPairFS::scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) {
  int* plusIndex = kmerPairIndices[PLUS];
  int* minusIndex = kmerPairIndices[MINUS];

  for (pos_t j=0; j<alignment->length; j++) {
    if (*plusIndex != -1)
      plusScores[j] += valueToWeight[*plusIndex];
    if (*minusIndex != -1)
      minusScores[j] += valueToWeight[*minusIndex];
    plusIndex++;
    minusIndex++;
  }
}
