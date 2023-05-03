#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " <target sequence file> <list of AXT files> <list of species names>" << endl;
    cerr << "List of species names should include target" << endl << endl;
    return 0;
  }

  if (argc % 2 != 1) {
    cerr << "Error, number of arguments must be even" << endl;
    return 0;
  }
  int numSequences = (argc - 1)/2;

  //read in list of sequences from command line
  //first one must be target sequence
  vector<string> sequences;
  for (int i=argc-numSequences; i<argc; i++) {
    sequences.push_back(argv[i]);
  }
  
  cerr << "Sequences in multiple alignment: ";
  for (int i=0; i<sequences.size(); i++)
    cerr << sequences[i] << " ";
  cerr << endl;

  //print the header
  cout << ">";
  for (int i=0; i<sequences.size(); i++) {
    cout << sequences[i];
    if (i != sequences.size()-1)
      cout << " ";
    else
      cout << endl;
  }

  // read in chromosome sequence
  cerr << "Reading chromosome sequence..." << endl;
  string line;
  string chr_seq = "";
  fstream chr_fs (argv[1], ios::in);
  getline (chr_fs, line);  /* skip header */
  while (getline (chr_fs, line)) {
    chr_seq += line;	
  }
  chr_fs.close();
  cerr << argv[1] << " length is " << chr_seq.length() << endl;

  //print the target sequence
  //mask any characters that are not bases
  for (int j=0; j<chr_seq.length(); j++) {
    if (chr_seq[j] != 'A' && chr_seq[j] != 'C' && chr_seq[j] != 'G' && chr_seq[j] != 'T' &&
	chr_seq[j] != 'a' && chr_seq[j] != 'c' && chr_seq[j] != 'g' && chr_seq[j] != 't')
      chr_seq[j] = 'N';
  }
  cerr << "Writing target sequence" << endl;
  cout << chr_seq << endl;
  
  // process alignments for each species
  for (int i=1; i<sequences.size(); i++) {
    cerr << "Processing alignment of " << sequences[i] << ": " << argv[i+1] << endl;

    char* alignmentRow = new char[chr_seq.length()];
   
    //initialize to all unaligned
    for (int j=0; j<chr_seq.length(); j++)
      alignmentRow[j] = '.';

    int blocks = 0;
    fstream axtStream (argv[i+1], ios::in);
    while (getline (axtStream, line)) {
      if (line == "" || line[0] == '#') continue;

      istringstream lineStream(line);
      int number;
      string targetChr;
      unsigned long targetStart;
      unsigned long targetEnd;
      string informantChr;
      unsigned long informantStart;
      unsigned long informantEnd;
      string strand;
      double score;
      lineStream >> number;
      lineStream >> targetChr;
      lineStream >> targetStart;
      lineStream >> targetEnd;
      lineStream >> informantChr;
      lineStream >> informantStart;
      lineStream >> informantEnd;
      lineStream >> strand;
      lineStream >> score;

      string targetSequence;
      string informantSequence;
      axtStream >> targetSequence;
      axtStream >> informantSequence;

      unsigned long targetPos = targetStart - 1;  // 0-based vs. 1-based coordinates
      for (int i=0; i<targetSequence.length(); i++) {
	if (targetSequence[i] != '-') {
	  alignmentRow[targetPos] = informantSequence[i];
	  targetPos++;
	}
      }

      blocks++;
      //if (blocks % 10000 == 0)
      //cerr << "Placed " << blocks << " blocks" << endl;
    }
    axtStream.close();

    //post-process and print out alignment row
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
