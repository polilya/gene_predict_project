#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <cassert>

#include "AlignmentBlock.h"

using namespace std;

int main(int argc, char** argv) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " <target sequence file> <MAF file> <sequence list>" << endl;
    return 0;
  }

  //read in list of sequences from command line
  //first one must be target sequence
  vector<string> sequences;
  for (int i=3; i<argc; i++) {
    sequences.push_back(argv[i]);
  }
  
  cerr << "Sequences in multiple alignment: ";
  for (int i=0; i<sequences.size(); i++)
    cerr << sequences[i] << " ";
  cerr << endl;

  string line;

  // read in chromosome sequence
  cerr << "Reading chromosome sequence..." << endl;
  string chr_seq = "";
  fstream chr_fs (argv[1], ios::in);
  getline (chr_fs, line);  /* skip header */
  while (getline (chr_fs, line)) {
    chr_seq += line;	
  }
  chr_fs.close();
  cerr << argv[1] << " length is " << chr_seq.length() << endl;
  

  //process MAF file
  vector<AlignmentBlock> blocks;
  int targetPosition;
  string targetSequence;
  std::string::const_iterator start, end;

  fstream maf_fs (argv[2], ios::in);
  while (getline (maf_fs, line)) {
    if (line[0] == 'a') {
      //start new alignment block
      if (blocks.size() != 0)
	blocks.back().compressBlock();
      blocks.push_back(AlignmentBlock());
      if (blocks.size() % 10000 == 0)
	cerr << "Reading block " << blocks.size() << endl;
      istringstream lineStream(line);
      string scoreString;
      lineStream >> scoreString;
      lineStream >> scoreString;
      int equalsPos = scoreString.find('=', 0);
      blocks.back().score = atof(scoreString.substr(equalsPos+1).c_str());
    }
    else if (line[0] == 's') {
      //alignment line
      istringstream lineStream(line);
      string s;
      string speciesString;
      int startPos;
      string sequence;
      lineStream >> s;
      lineStream >> speciesString;
      lineStream >> startPos;
      lineStream >> s;
      lineStream >> s;
      lineStream >> s;
      lineStream >> sequence;

      int dotPos = speciesString.find('.',0);
      string species = speciesString.substr(0, dotPos);

      /*
      cerr << "Read line.  Species = " << species << ", start = " 
	   << startPos << ", sequence = " << sequence << endl << endl;
      */

      if (species == sequences[0]) {
	blocks.back().start = startPos;
	blocks.back().length = sequence.length();
	blocks.back().addSequence(0, sequence);
      }
      else {
	//figure out what species this is
	for (int i=1; i<sequences.size(); i++) {
	  if (sequences[i] == species) {
	    //convert leading and trailing gaps to unaligned characters	 
	    for (int j=0; sequence[j] == '-' && j < sequence.length(); j++)
	      sequence[j] = '.';
	    for (int j=sequence.length()-1; sequence[j] == '-' && j>=0; j--)
	      sequence[j] = '.';
	    blocks.back().addSequence(i, sequence);
	  }
	}
      }
    }
  }
  maf_fs.close();

  //sort alignment blocks by score
  sort(blocks.begin(), blocks.end());

  //print the header
  cout << ">";
  for (int i=0; i<sequences.size(); i++) {
    cout << sequences[i];
    if (i != sequences.size()-1)
      cout << " ";
    else
      cout << endl;
  }

  //print the target sequence
  //mask any characters that are not bases
  for (int j=0; j<chr_seq.length(); j++) {
    if (chr_seq[j] != 'A' && chr_seq[j] != 'C' && chr_seq[j] != 'G' && chr_seq[j] != 'T' &&
	chr_seq[j] != 'a' && chr_seq[j] != 'c' && chr_seq[j] != 'g' && chr_seq[j] != 't')
      chr_seq[j] = 'N';
  }
  cerr << "Writing target sequence" << endl;
  cout << chr_seq << endl;

  //now create and print each row of the alignment
  for (int i=1; i<sequences.size(); i++) {
    char* alignmentRow = new char[chr_seq.length()+1];
    alignmentRow[chr_seq.length()] = '\0';
   
    //initialize to all unaligned
    for (int j=0; j<chr_seq.length(); j++)
      alignmentRow[j] = '.';

    //place blocks
    for (int j=0; j<blocks.size(); j++) {
      if (blocks[j].start + blocks[j].length > chr_seq.length()) {
	cerr << "Warning: Block starting at " << blocks[j].start << " has length " << blocks[j].length 
	     << ", which means it off the end of the chromosome (length " << chr_seq.length() << ")" << endl;
	continue;
      }

      int blockRow = -1;
      for (int k=0; k<blocks[j].species.size(); k++) {
	if (blocks[j].species[k] == i)
	  blockRow = k;
      }
      if (blockRow != -1) {
	unsigned long uncompressedLength = blocks[j].length;
	char* targetSequence = new char[blocks[j].length];
	char* alignedSequence = new char[blocks[j].length];
	uncompress((Bytef*)targetSequence, &uncompressedLength, 
		   (Bytef*)blocks[j].sequences[0], blocks[j].compressedSequenceLengths[0]);
	uncompress((Bytef*)alignedSequence, &uncompressedLength, 
		   (Bytef*)blocks[j].sequences[blockRow], blocks[j].compressedSequenceLengths[blockRow]);

	int index = blocks[j].start;
	for (int k=0; k<blocks[j].length; k++) {
	  if (targetSequence[k] != '-') {
	    alignmentRow[index] = alignedSequence[k];
	    index++;
	  }
	}

	delete[] targetSequence;
	delete[] alignedSequence;
      }
    }

    //post-process alignment row to convert all characters to uppercase
    //and change any unusual characters to N
    for (int j=0; j<chr_seq.length(); j++) {
      switch (alignmentRow[j]) {
      case 'a':
	alignmentRow[j] = 'A';
	break;
      case 'c':
	alignmentRow[j] = 'C';
	break;
      case 'g':
	alignmentRow[j] = 'G';
	break;
      case 't':
	alignmentRow[j] = 'T';
	break;
      case 'A':
	alignmentRow[j] = 'A';
	break;
      case 'C':
	alignmentRow[j] = 'C';
	break;
      case 'G':
	alignmentRow[j] = 'G';
	break;
      case 'T':
	alignmentRow[j] = 'T';
	break;
      case '-':
	alignmentRow[j] = '_';
	break;
      case '.':
	alignmentRow[j] = '.';
	break;
      default:
	alignmentRow[j] = 'N';
      }
    }

    cerr << "Writing alignment row " << i << endl;
    cout << alignmentRow << endl;
    delete[] alignmentRow;
  }

  return 0;
}
