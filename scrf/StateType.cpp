#include "StateType.h"
#include "DNAKmerFS.h"
#include "DNAKmerPairFS.h"

StateType::StateType() {
  lengthFeatures = NULL;
}

StateType::~StateType() {
}

void StateType::scorePositions() {
  if (sequenceFeatureSets.size() == 0) return;  /* nothing to do */

  pos_t length = sequenceFeatureSets[0]->alignment->length;

  plusScores.clear();
  minusScores.clear();
  plusScores.resize(length, 0);
  minusScores.resize(length, 0);
  for (int i=0; i<sequenceFeatureSets.size(); i++)
    sequenceFeatureSets[i]->scorePositions(plusScores, minusScores);
}

void StateType::freePositionScores() {
  vector<weight_t>().swap(plusScores);
  vector<weight_t>().swap(minusScores);
}

int StateType::findMatchingDNAKmerFS(int k, string parameterization, string seqName) {
  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    if (sequenceFeatureSets[i]->type == DNAKMER) {
      DNAKmerFS* dkfs = (DNAKmerFS*) sequenceFeatureSets[i];
      if (dkfs->k == k && dkfs->parameterization == parameterization && dkfs->seqName == seqName)
	return i;
    }
  }

  return -1;
}

int StateType::findMatchingDNAKmerPairFS(int k, string parameterization, string firstSeqName,
					 string secondSeqName) {
  for (int i=0; i<sequenceFeatureSets.size(); i++) {
    if (sequenceFeatureSets[i]->type == DNAKMERPAIR) {
      DNAKmerPairFS* dkpfs = (DNAKmerPairFS*) sequenceFeatureSets[i];
      if (dkpfs->k == k && dkpfs->parameterization == parameterization && 
	  dkpfs->firstSeqName == firstSeqName && dkpfs->secondSeqName == secondSeqName)
	return i;
    }
  }

  return -1;
}
