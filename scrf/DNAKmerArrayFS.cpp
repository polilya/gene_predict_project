#include <math.h>
#include "DNAKmerArrayFS.h"
#include "KmerIterator.h"

DNAKmerArrayFS::DNAKmerArrayFS(int k, int length, int offset, 
			       string parameterization, string seqName) {
  this->type = DNAKMERARRAY;
  this->k = k;
  this->length = length;
  this->offset = offset;
  this->parameterization = parameterization;
  this->seqName = seqName;

  int numValues = length * intpow(DNA_CHARS, k);
  valueToParam.resize(numValues);
  valueToWeight.resize(numValues);
  valueToGlobalParamID.resize(numValues);

  if (parameterization == "Standard") {
    int paramsPerPosition = intpow(4, k) + k + 1;
    int valuesPerPosition = intpow(DNA_CHARS, k);
    numParams = paramsPerPosition * length;
    
    KmerIterator kmerIter(k, "DNA");
    for (kmerIter.init(); ! kmerIter.finished; kmerIter.next()) {
      int numberOfNs = 0;
      bool containsGapOrUnaligned = false;
      for (int i=0; i<k; i++) {
	if (kmerIter.kmer[i] == 'N')
	  numberOfNs++;
	if (kmerIter.kmer[i] == '_' || kmerIter.kmer[i] == '.')
	  containsGapOrUnaligned = true;
      }
      int paramNum;
      if (containsGapOrUnaligned)
	paramNum = paramsPerPosition - 1;
      else if (numberOfNs > 0)
	paramNum = paramsPerPosition - 1 - numberOfNs;
      else {
	paramNum = 0;
	for (int i=0; i<k; i++)
	  paramNum += intpow(4, i) * kmerIter.counter[k - i - 1];
      }
      for (int i=0; i<length; i++) 
	valueToParam[kmerIter.index + i * valuesPerPosition] = paramNum + i * paramsPerPosition;
    }
  }
  else
    fatalError("Unknown parameterization for DNAKmerArrayFS");
}

void DNAKmerArrayFS::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  seqID = alignment->getSpeciesID(seqName);
  kmerIndices[PLUS] = kmerIndices[NONE] = alignment->kmerIndices[k][2 * seqID];
  kmerIndices[MINUS] = alignment->kmerIndices[k][2 * seqID + 1];
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;
}

void DNAKmerArrayFS::runningSumHelper(vector< vector<weight_t> >& runningSumPlus,
				      vector< vector<weight_t> >& runningSumMinus) {}

void DNAKmerArrayFS::scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) {}
