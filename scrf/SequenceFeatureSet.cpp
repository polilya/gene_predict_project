#include "SequenceFeatureSet.h"
#include "DNAKmerFS.h"
#include "DNAKmerPairFS.h"
#include "DNAKmerArrayFS.h"
#include "DNAKmerPairArrayFS.h"
#include "SVMFS.h"
#include "ESTPositionFS.h"
#include "ESTTransitionFS.h"
#include <set>

/* fills the denominators set with the global parameter indices of all
   parameters whose denominators should be incremented whenever value
   index is observed */
void SequenceFeatureSet::getDenominators(int index, set<int>& denominators) {
  denominators.clear();

  if (type == DNAKMER) {
    index = index - index % DNA_CHARS;
    for (int i=0; i<DNA_CHARS; i++) {
      denominators.insert(valueToGlobalParamID[index + i]);
    }
  }
  else if (type == DNAKMERPAIR) {
    index = index - index % DNA_CHARS;
    for (int i=0; i<DNA_CHARS; i++) {
      denominators.insert(valueToGlobalParamID[index + i]);
    }
  }
  else if (type == DNAKMERARRAY) {
    index = index - index % DNA_CHARS;
    for (int i=0; i<DNA_CHARS; i++) {
      denominators.insert(valueToGlobalParamID[index + i]);
    }
  }
  else if (type == DNAKMERPAIRARRAY) {
    index = index - index % DNA_CHARS;
    for (int i=0; i<DNA_CHARS; i++) {
      denominators.insert(valueToGlobalParamID[index + i]);
    }
  }
}

void SequenceFeatureSet::positionInitialGuessHelper(pos_t pos, strand_t strand,
						    vector<weight_t>& denominators) {
  int index;

  if (type == MASKING) {
    denominators[valueToGlobalParamID[0]] += 1.0;
  }
  else if (type == DNAKMER) {
    DNAKmerFS* dkfs = (DNAKmerFS*)this;
    
    if (strand == PLUS || strand == NONE) {
      if (pos > alignment->length - dkfs->k) return; 
      index = alignment->getDNAKmerIndex(2 * dkfs->seqID, pos, dkfs->k);
    }
    else if (strand == MINUS) {
      if (alignment->length - pos - 1 > alignment->length - dkfs->k) return;
      index = alignment->getDNAKmerIndex(2 * dkfs->seqID + 1, alignment->length - pos - 1, dkfs->k);
    }
    set<int> possibleParams;
    index = index - index % DNA_CHARS;
    for (int i=0; i<DNA_CHARS; i++) {
      possibleParams.insert(valueToGlobalParamID[index + i]);
    }
    for (set<int>::iterator iter = possibleParams.begin(); iter != possibleParams.end(); ++iter) {
      denominators[*iter] += (weight_t)1;
    }
  }
  else if (type == DNAKMERPAIR) {
    DNAKmerPairFS* dkpfs = (DNAKmerPairFS*)this;
    
    if (strand == PLUS || strand == NONE) {
      if (pos > alignment->length - dkpfs->k) return; 
      index = alignment->getDNAKmerPairIndex(2 * dkpfs->firstSeqID, 2 * dkpfs->secondSeqID, pos, dkpfs->k);
    }
    else if (strand == MINUS) {
      if (alignment->length - pos - 1 > alignment->length - dkpfs->k) return;
      index = alignment->getDNAKmerPairIndex(2 * dkpfs->firstSeqID + 1, 2 * dkpfs->secondSeqID + 1, 
				  alignment->length - pos - 1, dkpfs->k);
    }
    set<int> possibleParams;
    index = index - index % DNA_CHARS;
    for (int i=0; i<DNA_CHARS; i++) {
      possibleParams.insert(valueToGlobalParamID[index + i]);
    }
    for (set<int>::iterator iter = possibleParams.begin(); iter != possibleParams.end(); ++iter) {
      denominators[*iter] += (weight_t)1;
    }
  }
  else if (type == DNAKMERARRAY) {
    DNAKmerArrayFS* dkafs = (DNAKmerArrayFS*)this;
    pos_t start;
    int row;

    if (strand == PLUS || strand == NONE) {
      start = pos - dkafs->offset;
      row = 2 * dkafs->seqID;
    }
    else if (strand == MINUS) {
      start = alignment->length - 1 - pos - dkafs->offset + 1;
      row = 2 * dkafs->seqID + 1;
    }

    if (start >= 0 && start + dkafs->length + dkafs->k <= alignment->length) {
      for (int i=0; i<dkafs->length; i++) {
	int index = alignment->getDNAKmerIndex(row, start + i, dkafs->k);
	index = index - index % DNA_CHARS;
	for (int j=0; j<DNA_CHARS; j++) {
	  denominators[valueToGlobalParamID[index + j + i * intpow(DNA_CHARS, dkafs->k)]] += (weight_t) 1;
	}
      }
    }
  }
  else if (type == DNAKMERPAIRARRAY) {
    DNAKmerPairArrayFS* dkpafs = (DNAKmerPairArrayFS*)this;
    pos_t start;
    int firstRow;
    int secondRow;

    if (strand == PLUS || strand == NONE) {
      start = pos - dkpafs->offset;
      firstRow = dkpafs->firstSeqID * 2;
      secondRow = dkpafs->secondSeqID * 2;
    }
    else if (strand == MINUS) {
      start = alignment->length - 1 - pos - dkpafs->offset + 1;
      firstRow = dkpafs->firstSeqID * 2 + 1;
      secondRow = dkpafs->secondSeqID * 2 + 1;
    }

    if (start >= 0 && start + dkpafs->length + dkpafs->k <= alignment->length) {
      for (int i=0; i<dkpafs->length; i++) {
	int index = alignment->getDNAKmerPairIndex(firstRow, secondRow, start + i, dkpafs->k) +
	  i*intpow(DNA_CHARS, 2*dkpafs->k);
	set<int> possibleParams;
	index = index - index % DNA_CHARS;
	for (int j=0; j<DNA_CHARS; j++) {
	  possibleParams.insert(valueToGlobalParamID[index + j]);
	}
	for (set<int>::iterator iter = possibleParams.begin(); iter != possibleParams.end(); ++iter) {
	  denominators[*iter] += (weight_t)1;
	}
      }
    }
  }
  else if (type == SVM) {

  }
  else if (type == ESTPOSITION) {
    for (int i=0; i<EST_CHARS; i++) {
      denominators[valueToGlobalParamID[i]] += 1.0;
    }
  }
  else if (type == ESTTRANSITION) {
    for (int i=0; i<EST_CHARS; i++) {
      denominators[valueToGlobalParamID[i]] += 1.0;
    }
  }
  else
    fatalError("positionInitialGuessHelper called with unknown SFS type");
}

void SequenceFeatureSet::positionExpectationHelper(pos_t pos, strand_t strand, weight_t prob, 
						   vector<weight_t>& counts) {
  int index;
  
  if (type == MASKING) {
    switch (alignment->sequenceArray[0][pos]) {
    case 'a':
    case 'c':
    case 'g':
    case 't':
      counts[valueToGlobalParamID[0]] += prob;
    }
  }
  else if (type == DNAKMER) {
    DNAKmerFS* dkfs = (DNAKmerFS*)this;
    index = dkfs->kmerIndices[strand][pos];
    if (index != -1)
      counts[valueToGlobalParamID[index]] += prob;   
  }
  else if (type == DNAKMERPAIR) {
    DNAKmerPairFS* dkpfs = (DNAKmerPairFS*)this;    
    index = dkpfs->kmerPairIndices[strand][pos];
    if (index != -1)
      counts[valueToGlobalParamID[index]] += prob;
  }
  else if (type == DNAKMERARRAY) {
    DNAKmerArrayFS* dkafs = (DNAKmerArrayFS*)this;

    pos_t start = strand == MINUS ? alignment->length - 1 - pos - dkafs->offset + 1 : pos - dkafs->offset;
    if (start >= 0 && start + dkafs->length + dkafs->k <= alignment->length) {
      for (int i=0; i<dkafs->length; i++) {
	index = dkafs->kmerIndices[strand][start + i];
	index += i * DNA_INDEX_COEFF[dkafs->k];
	counts[valueToGlobalParamID[index]] += prob;
      }
    }
  }
  else if (type == DNAKMERPAIRARRAY) {
    DNAKmerPairArrayFS* dkpafs = (DNAKmerPairArrayFS*)this;

    pos_t start = strand == MINUS ? alignment->length - 1 - pos - dkpafs->offset + 1 : pos - dkpafs->offset;
    if (start >= 0 && start + dkpafs->length + dkpafs->k <= alignment->length) {
      for (int i=0; i<dkpafs->length; i++) {
	index = dkpafs->kmerPairIndices[strand][start + i];
	index += i * DNA_INDEX_COEFF[2 * dkpafs->k];
	counts[valueToGlobalParamID[index]] += prob;
      }
    }
  }
  else if (type == SVM) {
    SVMFS* svmfs = (SVMFS*)this;

    int bin = svmfs->bins[strand][pos];
    float coeff = svmfs->coeffs[strand][pos];
    counts[valueToGlobalParamID[bin]] += prob * coeff;
    counts[valueToGlobalParamID[bin+1]] += prob * (1 - coeff);
  }
  else if (type == ESTPOSITION) {
    index = EST_TO_INDEX[estSeq->sequence[pos]];
    counts[valueToGlobalParamID[index]] += prob;   
  }
  else if (type == ESTTRANSITION) {
    if (strand == PLUS || strand == NONE)
      index = EST_CHARS * EST_TO_INDEX[estSeq->sequence[pos-1]] + EST_TO_INDEX[estSeq->sequence[pos]];
    else
      index = EST_TO_INDEX[estSeq->sequence[pos-1]] + EST_CHARS * EST_TO_INDEX[estSeq->sequence[pos]];
    counts[valueToGlobalParamID[index]] += prob;
  }
  else
    fatalError("positionExpectationHelper called with unknown SFS type");
}
