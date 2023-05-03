#include "DNAKmerPairArrayFS.h"
#include "KmerIterator.h"

DNAKmerPairArrayFS::DNAKmerPairArrayFS(int k, int length, int offset,
				       string parameterization, string firstSeqName,
				       string secondSeqName) {
  this->type = DNAKMERPAIRARRAY;
  this->k = k;
  this->length = length;
  this->offset = offset;
  this->parameterization = parameterization;
  this->firstSeqName = firstSeqName;
  this->secondSeqName = secondSeqName;

  int numValues = length * intpow(DNA_CHARS, 2*k);
  valueToParam.resize(numValues);
  valueToWeight.resize(numValues);
  valueToGlobalParamID.resize(numValues);

  if (parameterization == "Standard") {
    int paramsPerPosition = intpow(4,2*k) + (intpow(2,k)-1) + 3 + 1;
    int valuesPerPosition = intpow(DNA_CHARS, 2*k);
    numParams = paramsPerPosition * length;

    KmerIterator kmerIter(2*k, "DNA");
    for (kmerIter.init(); ! kmerIter.finished; kmerIter.next()) {
      int firstKmerParam = kmerToParameter(kmerIter.kmer, k);
      int secondKmerParam = kmerToParameter(kmerIter.kmer + k, k);
      int paramNum;

      if (firstKmerParam >= intpow(4,k))
	paramNum = paramsPerPosition - 1;
      else if (secondKmerParam < intpow(4,k))
	paramNum = intpow(4,k) * firstKmerParam + secondKmerParam;
      else
	paramNum = intpow(4,2*k) + (secondKmerParam - intpow(4,k));
     
      for (int i=0; i<length; i++)
	valueToParam[kmerIter.index + i * valuesPerPosition] = paramNum + i * paramsPerPosition;
    }
  }
  else if (parameterization == "General") {
    numParams = numValues;
    for (int i=0; i<numValues; i++) {
      valueToParam[i] = i;
    }
  }
  else
    fatalError("Unknown parameterization for DNAKmerPairArrayFS");
}

//gets the parameter number of an unpaired k-mer
int DNAKmerPairArrayFS::kmerToParameter(char* kmer, int k) {
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

void DNAKmerPairArrayFS::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  firstSeqID = alignment->getSpeciesID(firstSeqName);
  secondSeqID = alignment->getSpeciesID(secondSeqName);
  kmerPairIndices[PLUS] = kmerPairIndices[NONE] = alignment->kmerPairIndices[k][2 * firstSeqID][2 * secondSeqID];
  kmerPairIndices[MINUS] = alignment->kmerPairIndices[k][2 * firstSeqID + 1][2 * secondSeqID + 1];
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;
}

void DNAKmerPairArrayFS::runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				      vector< vector<weight_t> >& runningSumMinus) {}

void DNAKmerPairArrayFS::scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) {}
