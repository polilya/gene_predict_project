#include <math.h>
#include "DNAKmerFS.h"
#include "KmerIterator.h"


DNAKmerFS::DNAKmerFS(int k, int frame, string parameterization, string seqName) {
  this->type = DNAKMER;
  this->k = k;
  this->frame = frame;
  this->parameterization = parameterization;
  this->seqName = seqName;
  seqID = -1;  //unknown as yet

  int numValues = intpow(DNA_CHARS, k);
  valueToParam.resize(numValues);
  valueToWeight.resize(numValues);
  valueToGlobalParamID.resize(numValues);

  if (parameterization == "Standard") {
    numParams = intpow(4, k) + k + 1;

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
	paramNum = numParams - 1;
      else if (numberOfNs > 0)
	paramNum = numParams - 1 - numberOfNs;
      else {
	paramNum = 0;
	for (int i=0; i<k; i++)
	  paramNum += intpow(4, i) * kmerIter.counter[k - i - 1];
      }
      valueToParam[kmerIter.index] = paramNum;
    }
  }
  else {
    cerr << "Parameterization " << parameterization << " not recognized." << endl;
    fatalError("Unknown parameterization for DNAKmerFS");
  }
}

void DNAKmerFS::setSequences(AlignmentSequence* alignment, NumericSequence* numericSeq, ESTSequence* estSeq) {
  this->alignment = alignment;
  seqID = alignment->getSpeciesID(seqName);
  kmerIndices[PLUS] = kmerIndices[NONE] = alignment->kmerIndices[k][2 * seqID];
  kmerIndices[MINUS] = alignment->reverseKmerIndices[k][2 * seqID + 1];
  this->numericSeq = numericSeq;
  this->estSeq = estSeq;
}

void DNAKmerFS::scorePositions(vector<weight_t>& plusScores, vector<weight_t>& minusScores) {
  int* plusIndex = kmerIndices[PLUS];
  int* minusIndex = kmerIndices[MINUS];

  for (pos_t j=0; j<alignment->length; j++) {
    if (*plusIndex != -1)
      plusScores[j] += valueToWeight[*plusIndex];
    if (*minusIndex != -1)
      minusScores[j] += valueToWeight[*minusIndex];
    plusIndex++;
    minusIndex++;
  }
}
