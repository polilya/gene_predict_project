#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " <sequence file> <RepeatMasker .out file>" << endl;
    return 0;
  }

  string header;
  string line;

  // read in chromosome sequence
  cerr << "Reading chromosome sequence..." << endl;
  string chr_seq = "";
  fstream chr_fs (argv[1], ios::in);
  getline (chr_fs, header);  // header
  while (getline (chr_fs, line)) {
    chr_seq += line;	
  }
  chr_fs.close();
  cerr << argv[1] << " length is " << chr_seq.length() << endl;
  
  //convert all lowercase a,c,g,t to uppercase and any non-ACGT characters to N
  for (int j=0; j<chr_seq.length(); j++) {
    switch (chr_seq[j]) {
    case 'a':
      chr_seq[j] = 'A';
      break;
    case 'c':
      chr_seq[j] = 'C';
      break;
    case 'g':
      chr_seq[j] = 'G';
      break;
    case 't':
      chr_seq[j] = 'T';
      break;
    case 'A':
      chr_seq[j] = 'A';
      break;
    case 'C':
      chr_seq[j] = 'C';
      break;
    case 'G':
      chr_seq[j] = 'G';
      break;
    case 'T':
      chr_seq[j] = 'T';
      break;
    default:
      chr_seq[j] = 'N';
    }
  }
    
  //process RepeatMasker .out file
  fstream rm_fs (argv[2], ios::in);

  //skip headers
  getline(rm_fs, line);
  getline(rm_fs, line);
  getline(rm_fs, line);

  while (getline (rm_fs, line)) {
    //cerr << "Line: " << line << endl;
    int begin;
    int end;
    string family;

    istringstream lineStream (line, istringstream::in);
    string s;
    lineStream >> s;  //SW score
    lineStream >> s;  //perc div.
    lineStream >> s;  //perc del.
    lineStream >> s;  //perc ins.
    lineStream >> s;  //query sequence

    lineStream >> begin;
    begin -= 1;  //adjust for 1-based coordinates

    lineStream >> end;
    end -= 1;  //adjust for 1-based coordinates

    lineStream >> s;  //(left)
    lineStream >> s;  //strand
    lineStream >> s;  //matching repeat

    lineStream >> family;

    if (family != "Simple_repeat" && family != "Low_complexity") {
      //change all ACGT characters to lowercase
      //cerr << "Masking repeat of type " << family << " from " << begin << " to " << end << endl;
      for (int j=begin; j<=end; j++) {
	switch(chr_seq[j]) {
	case 'A':
	  chr_seq[j] = 'a';
	  break;
	case 'C':
	  chr_seq[j] = 'c';
	  break;
	case 'G':
	  chr_seq[j] = 'g';
	  break;
	case 'T':
	  chr_seq[j] = 't';
	  break;
	}
      }
    }
  }
 
  //print the header
  cout << header << endl;

  //print the masked sequence
  cout << chr_seq << endl;

  return 0;
}
