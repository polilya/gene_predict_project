#include <fstream>
#include <sstream>
#include <zlib.h>
#include "AlignmentSequence.h"
#include "Globals.h"

#include "MemoryDebug.h"

/* Reads in a .align file */
AlignmentSequence::AlignmentSequence(string alignFilename, set<string>& speciesToRead) {
  const int BUFFER_SIZE = 1048576;  //2^20

  string alignment = "";
  char buffer[BUFFER_SIZE];

  gzFile alignFile;
  alignFile = gzopen (alignFilename.c_str(), "rb");
  if (alignFile == NULL) {
    cerr << "Could not open " << alignFilename << endl;
    fatalError("Could not open alignment file");
  }
  while (! gzeof(alignFile)) {
    int bytesRead = gzread (alignFile, buffer, BUFFER_SIZE-1);
    buffer[bytesRead] = '\0';
    alignment += buffer;
  }
  gzclose (alignFile);

  istringstream alignStream(alignment);
 
  string s;

  /* parse the header */
  char headerArray[STRLEN];
  alignStream.getline(headerArray, 1024);
  string header(headerArray);
  header = header.substr(1, s.length() - 1);  /* trim off the leading '>' character */
  istringstream headerStream(header);
  while (! headerStream.eof()) {
    headerStream >> s;
    //if (! headerStream.eof())
      species.push_back(s);	
  }

  numberOfSequences = 2 * (species.size() + 1);
  
  sequenceArray = new char*[numberOfSequences];
  for (int i = 0; i < numberOfSequences; i++)
    sequenceArray[i] = NULL;
  
  compressedSequenceArray = new char*[species.size()];
  for (int i = 0; i < species.size(); i++)
    compressedSequenceArray[i] = NULL;
  compressedSequenceLengths.resize(species.size());

  /* read in the sequences */
  int i = 0;
  while (! alignStream.eof()) {
    alignStream >> s;
    if (! alignStream.eof()) {
      length = s.length();
      if (speciesToRead.find(species[i]) != speciesToRead.end()) {
	sequenceArray[2 * i] = new char[length];
	memcpy(sequenceArray[2*i], s.c_str(), length);
      }
      i++;
    }
  }
  if (i < species.size()) {
    cerr << alignFilename << " contained fewer lines than indicated by its header" << endl;
    fatalError("Alignment file read error");
  }

  compressSequences();
  freeUncompressedSequences();

  kmerIndices = NULL;
  kmerPairIndices = NULL;
  reverseKmerIndices = NULL;
  reverseKmerPairIndices = NULL;
}

AlignmentSequence::~AlignmentSequence() {
  freeKmerIndices();
  freeUncompressedSequences();
  freeCompressedSequences();
  delete[] sequenceArray;
  delete[] compressedSequenceArray;
}

void AlignmentSequence::write(string outFilename) {
  ofstream outStream(outFilename.c_str());
  outStream << ">";
  for (int i=0; i<species.size(); i++) {
    outStream << species[i];
    if (i == species.size() - 1)
      outStream << endl;
    else
      outStream << " ";
  }
  for (int i=0; i<species.size(); i++)
    outStream << sequenceArray[2 * i] << endl;
}

int AlignmentSequence::getSpeciesID(string s) {
  for (int i=0; i<species.size(); i++) {
    if (species[i] == s)
      return i;
  }
  cerr << "No match for species " << s << endl;
  fatalError("Unable to find species in getSpeciesID");
  return -1;
}

void AlignmentSequence::reverseComplement(char* dna, char* rc) {
  for (pos_t i=0; i<length; i++)
    rc[length - i - 1] = INDEX_TO_DNA[DNA_INDEX_COMPLEMENT[DNA_TO_INDEX[dna[i]]]];
}

void AlignmentSequence::createReverseComplements() {
  for (int i = 0; i < species.size(); i++) {
    if (sequenceArray[2*i] != NULL) {
      sequenceArray[2*i + 1] = new char[length];
      reverseComplement(sequenceArray[2*i], sequenceArray[2*i + 1]);
    }
  }
}

void AlignmentSequence::freeReverseComplements() {
  for (int i = 0; i < species.size(); i++) {
    if (sequenceArray[2*i+1] != NULL) {
      delete[] sequenceArray[2*i + 1];
      sequenceArray[2*i+1] = NULL;
    }
  }

}

//converts all lowercase bases in the target sequence to N
void AlignmentSequence::maskLowercase() {
  for (pos_t i=0; i<length; i++) {
    if (sequenceArray[0][i] == 'a' || sequenceArray[0][i] == 'c' || 
	sequenceArray[0][i] == 'g' || sequenceArray[0][i] == 't')
      sequenceArray[0][i] = 'N';
  }
}


void AlignmentSequence::compressSequences() {
  freeCompressedSequences();

  unsigned long bufferSize = 2 * length;
  char* buffer = new char[bufferSize];

  for (int i=0; i<species.size(); i++) {
    if (sequenceArray[2*i] != NULL) {
      compressedSequenceLengths[i] = bufferSize;
      compress2((Bytef*)buffer, &(compressedSequenceLengths[i]), (Bytef*)sequenceArray[2*i], length, 1);
      compressedSequenceArray[i] = new char[compressedSequenceLengths[i]];
      memcpy(compressedSequenceArray[i], buffer, compressedSequenceLengths[i]);
    }
  }

  delete[] buffer;
}

void AlignmentSequence::uncompressSequences() {
  freeUncompressedSequences();

  unsigned long uncompressedLength = length;

  for (int i = 0; i < species.size(); i++) {
    if (compressedSequenceArray[i] != NULL) {
      sequenceArray[2*i] = new char[length];
      uncompress((Bytef*)sequenceArray[2*i], &uncompressedLength, 
		 (Bytef*)compressedSequenceArray[i], compressedSequenceLengths[i]);
    }
  }
}

void AlignmentSequence::freeUncompressedSequences() {
  for (int i=0; i<numberOfSequences; i++) {
    if (sequenceArray[i] != NULL) {
      delete[] sequenceArray[i];
      sequenceArray[i] = NULL;
    }
  }
}

void AlignmentSequence::freeCompressedSequences() {
  for (int i=0; i<species.size(); i++) {
    if (compressedSequenceArray[i] != NULL) {
      delete[] compressedSequenceArray[i];
      compressedSequenceArray[i] = NULL;
    }
  }
}


void AlignmentSequence::freeKmerIndices() {
  if (kmerIndices != NULL) {
    for (int i=0; i<=maxK; i++) {
      if (kmerIndices[i] != NULL) {
	for (int j=0; j<=2*maxSeq+1; j++) {
	  if (kmerIndices[i][j] != NULL) {
	    delete[] kmerIndices[i][j];
	    delete[] reverseKmerIndices[i][j];
	  }
	}
	delete[] kmerIndices[i];
	delete[] reverseKmerIndices[i];
      }
    }
    delete[] kmerIndices;
    delete[] reverseKmerIndices;
    kmerIndices = NULL;
    reverseKmerIndices = NULL;
  }

  if (kmerPairIndices != NULL) {
    for (int i=0; i<=maxPairK; i++) {
      if (kmerPairIndices[i] != NULL) {
	for (int j=0; j<=2*maxPairSeq+1; j++) {
	  if (kmerPairIndices[i][j] != NULL) {
	    for (int k=0; k<=2*maxPairSeq+1; k++) {
	      if (kmerPairIndices[i][j][k] != NULL) {
		delete[] kmerPairIndices[i][j][k];
		delete[] reverseKmerPairIndices[i][j][k];
	      }
	    }
	    delete[] kmerPairIndices[i][j];
	    delete[] reverseKmerPairIndices[i][j];
	  }
	}
	delete[] kmerPairIndices[i];
	delete[] reverseKmerPairIndices[i];
      }
    }
    delete[] kmerPairIndices;
    delete[] reverseKmerPairIndices;
    kmerPairIndices = NULL;
    reverseKmerPairIndices = NULL;
  }
}
