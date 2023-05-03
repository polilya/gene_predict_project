#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

using namespace std;

const int MAX_GAP = 1;
const float MIN_IDENTITY = 0.98;
const int MIN_INTRON_LENGTH = 20;

void printStats(char* estseq, string chr_seq) {
  char chars[] = {'U', 'S', 'I', 'N', 'C'};
  int charCounts[5];
  for (int i=0; i<5; i++)
    charCounts[i] = 0;
  for (int i=0; i<chr_seq.length(); i++) {
    switch (estseq[i]) {
    case 'U':
      charCounts[0]++;
      break;
    case 'S':
      charCounts[1]++;
      break;
    case 'I':
      charCounts[2]++;
      break;
    case 'N':
      charCounts[3]++;
      break;
    case 'C':
      charCounts[4]++;
      break;
    }
  }
  for (int i=0; i<5; i++)
    cerr << chars[i] << ": " << (double)charCounts[i]/chr_seq.length() << endl;
}


int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " <chromosome FASTA file> <UCSC est alignment file>" << endl;
    exit(0);
  }

  string s, line;
  int integer;

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

  char* estseq = new char[chr_seq.length()];
  for (int i=0; i<chr_seq.length(); i++)
    estseq[i] = 'N';

  //process EST alignments
  int estsProcessed = 0;
  fstream est_fs (argv[2], ios::in);
  while (getline (est_fs, line)) {
    vector<string> lineArray;
    istringstream lineStream(line);
    while (! lineStream.eof()) {
      lineStream >> s;
      lineArray.push_back(s);
    }

    int matches = atoi(lineArray[1].c_str());
    int length = atoi(lineArray[11].c_str());
    float identity = (float)matches / length;
    if (identity < MIN_IDENTITY) continue;

    vector<int> blatStarts;
    vector<int> blatSizes;
    replace(lineArray[19].begin(), lineArray[19].end(), ',', ' ');
    replace(lineArray[21].begin(), lineArray[21].end(), ',', ' ');
    istringstream blatStartStream(lineArray[21]);
    while (! blatStartStream.eof()) {
      blatStartStream >> integer;
      if (! blatStartStream.eof())
	blatStarts.push_back(integer);
    }
    istringstream blatSizeStream(lineArray[19]);
    while (! blatSizeStream.eof()) {
      blatSizeStream >> integer;
      if (! blatSizeStream.eof())
	blatSizes.push_back(integer);
    }
    assert(blatStarts.size() == blatSizes.size());
    assert(blatStarts.size() == atoi(lineArray[18].c_str()));

    //merge blocks separated by small gaps
    vector<int> blockStarts;
    vector<int> blockSizes;
    int i = 0;
    while (i < blatStarts.size()) {
      blockStarts.push_back(blatStarts[i]);
      int end = blatStarts[i] + blatSizes[i] - 1;
      for (int j=i+1; blatStarts[j] - end <= MAX_GAP+1 && j<blatStarts.size(); j++) {
	end = blatStarts[j] + blatSizes[j] - 1;
	i++;
      }
      blockSizes.push_back(end - blockStarts.back() + 1);
      i++;
    }

    /*
    cerr << "BLAT blocks" << endl;
    for (int i=0; i<blatStarts.size(); i++)
      cerr << blatStarts[i] << "\t" << (blatStarts[i] + blatSizes[i] - 1) << endl;
    cerr << "Blocks" << endl;
    for (int i=0; i<blockStarts.size(); i++)
      cerr << blockStarts[i] << "\t" << (blockStarts[i] + blockSizes[i] - 1) << endl;
    cerr << endl;
    */

    //check to make sure than any remaining gaps in the EST alignment look like plausible introns
    string strand = lineArray[9];
    bool badGap = false;
    for (int i=0; i<blockStarts.size()-1; i++) {
      int intronStart = blockStarts[i] + blockSizes[i];
      int intronEnd = blockStarts[i+1] - 1;
      int intronLength = intronEnd - intronStart + 1;
      if (intronLength < MIN_INTRON_LENGTH)
	badGap = true;
      if (strand == "+") {
	if (chr_seq[intronStart] != 'G' || (chr_seq[intronStart+1] != 'T' && chr_seq[intronStart+1] != 'C') ||
	    chr_seq[intronEnd-1] != 'A' || chr_seq[intronEnd] != 'G')
	  badGap = true;
      }
      else if (strand == "-") {
	if (chr_seq[intronStart] != 'C' || chr_seq[intronStart+1] != 'T' ||
	    (chr_seq[intronEnd-1] != 'A' && chr_seq[intronEnd-1] != 'G') || chr_seq[intronEnd] != 'C')
	  badGap = true;
      }
      else {
	cerr << "Bad strand for EST" << endl;
	exit(0);
      }
    }
    if (badGap) continue;

    bool splicedEST = blockStarts.size() > 1;
    for (int i=0; i<blockStarts.size(); i++) {
      //place exon block
      for (int j=blockStarts[i]; j<blockStarts[i] + blockSizes[i]; j++) {
	if (splicedEST) {
	  if (estseq[j] == 'N' || estseq[j] == 'S' || estseq[j] == 'U')
	    estseq[j] = 'S';
	  else
	    estseq[j] = 'C';
	}
	else {
	  if (estseq[j] == 'N')
	    estseq[j] = 'U';
	  //if unspliced EST conflicts with spliced EST evidence, just ignore it
	}
      }
      if (i != blockStarts.size() - 1) {
	//place intron block
	for (int j=blockStarts[i]+blockSizes[i]; j<blockStarts[i+1]; j++) {
	  if (estseq[j] == 'N' || estseq[j] == 'I' || estseq[j] == 'U')
	    estseq[j] = 'I';  //intron from spliced EST overrides unspliced EST evidence
	  else
	    estseq[j] = 'C';
	}
      }
    }
    estsProcessed++;
    if (estsProcessed % 10000 == 0) {
      cerr << "Processed " << estsProcessed << " ESTs" << endl;
    }
  }
  est_fs.close();

  printStats(estseq, chr_seq);

  cout << ">EST Sequence" << endl;
  cout << estseq << endl;

  return 0;
}

