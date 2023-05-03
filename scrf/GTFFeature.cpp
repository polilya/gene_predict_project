#include "GTFFeature.h"
#include <sstream>


void GTFFeature::read(ifstream& gtfStream) {
  string strandString;
  string scoreString;
  string fohString;

  gtfStream >> seqName >> source >> featureName >> start >> end >>
    scoreString >> strandString >> fohString;
  
  if (strandString == "+")
    strand = PLUS;
  else if (strandString == "-")
    strand = MINUS;
  else if (! gtfStream.eof()) {
    fatalError("Illegal strand encountered in GTF file");
  }

  /* score and foh are allowed to be undefined */
  if (scoreString == ".") {
    score = (weight_t)0;
    scoreDefined = false;
  }
  else {
    score = atof(scoreString.c_str());
    scoreDefined = true;
  }

  if (fohString == "0" || fohString == "1" || fohString == "2")
    foh = atoi(fohString.c_str());
  else
    foh = -1;

  
  /* read attributes */
  char attributeLine[STRLEN];
  gtfStream.getline(attributeLine, STRLEN, '\n');
  istringstream attributeLineStream(attributeLine);

  while (! attributeLineStream.eof()) {
    char attribute[STRLEN];
    string tag;
    string value;

    attributeLineStream.getline(attribute, STRLEN, ';');
    istringstream attributeStream(attribute);
    attributeStream >> tag;
    attributeStream >> value;
    
    /* delete the quotes from the tag */
    while (value.find('"', 0) != string::npos) {
      value.erase(value.find('"', 0), 1);
    }

    if (tag == "gene_id")
      geneID = value;
    else if (tag == "transcript_id")
      transcriptID = value;
  }
}

void GTFFeature::print(ostream& outStream) {
  outStream << seqName << "\t" << source << "\t" << featureName;
  //if (featureName == "CDS") outStream << "\t";
  outStream << "\t" << start << "\t" << end << "\t";

  if (scoreDefined) {
    ios_base::fmtflags oldFlags = outStream.flags();
    streamsize oldPrecision = outStream.precision(4);
    outStream << fixed;
    outStream << score;
    outStream.precision(oldPrecision);
    outStream.flags(oldFlags);
  }
  else
    outStream << ".";
  
  outStream << "\t";

  if (strand == PLUS)
    outStream << "+";
  else if (strand == MINUS)
    outStream << "-";
  else
    fatalError("Illegal strand in GTF feature");

  outStream << "\t";
  
  if (foh == -1)
    outStream << ".";
  else
    outStream << foh;

  outStream << "\t" << "gene_id \""
	    << geneID << "\"; " << "transcript_id \""
	    << transcriptID << "\";" << endl;   
}

bool GTFFeature::operator==(const GTFFeature& other) {
  return featureName == other.featureName && start == other.start && end == other.end;
}

bool GTFFeature::operator!=(const GTFFeature& other) {
  return !operator==(other);
}

bool GTFFeature::operator<(const GTFFeature& other) const {
  if (start != other.start)
    return start < other.start;
  else if (featureName != other.featureName) {
    //sort by feature type
    int rank = 0;
    int otherRank = 0;

    if (featureName == "5UTR")
      rank = 0;
    else if (featureName == "start_codon")
      rank = 1;
    else if (featureName ==  "CDS")
      rank = 2;
    else if (featureName == "stop_codon")
      rank = 3;
    else if (featureName == "3UTR")
      rank = 4;
    else
      fatalError("Unknown feature name in GTF file");

    if (other.featureName == "5UTR")
      otherRank = 0;
    else if (other.featureName == "start_codon")
      otherRank = 1;
    else if (other.featureName ==  "CDS")
      otherRank = 2;
    else if (other.featureName == "stop_codon")
      otherRank = 3;
    else if (other.featureName == "3UTR")
      otherRank = 4;
    else
      fatalError("Unknown feature name in GTF file");

    if (rank != otherRank) {
      if (strand == PLUS)
	return rank < otherRank;
      else if (strand == MINUS)
	return rank > otherRank;
    }
  }
  else if (end != other.end)
    return end < other.end;
  
  return false;  /* tie */
}
