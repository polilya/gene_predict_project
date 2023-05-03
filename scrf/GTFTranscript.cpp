#include "GTFTranscript.h"

void GTFTranscript::setAttributesFromFeatures() {
  id = features.front().transcriptID;
  strand = features.front().strand;
  start = features.front().start;
  end = features.back().end;
}

void GTFTranscript::print(ostream& outStream) {
  for (int i=0; i<features.size(); i++)
    features[i].print(outStream);
}

bool GTFTranscript::overlaps(GTFTranscript& other) {
  if (start >= other.start && start <= other.end)
    return true;
  else if (end >= other.start && end <= other.end)
    return true;
  else if (start <= other.start && end >= other.end)
    return true;
  else
    return false;
}


void GTFTranscript::setFohs() {
  int foh, toh = 0;
  if (strand == PLUS) {
    for (int i=0; i<features.size(); i++) {
      if (features[i].featureName == "CDS") {
	int length = features[i].end - features[i].start + 1;
	foh = (3 - toh) % 3;
	toh = TOH_LOOKUP[foh][length % 3];
	features[i].foh = foh;
      }
    }
  }
  else {
    for (int i=features.size()-1; i >= 0; i--) {
      if (features[i].featureName == "CDS") {
	int length = features[i].end - features[i].start + 1;
	foh = (3 - toh) % 3;
	toh = TOH_LOOKUP[foh][length % 3];
	features[i].foh = foh;
      }
    }
  }
}

void GTFTranscript::setExonTypes() {
  //first sort features
  sort(features.begin(), features.end());

  for (int i=0; i<features.size(); i++) {
    if (features[i].featureName == "5UTR") {
      if (strand == PLUS) {
	if (i == 0) {	  
	  if (i+1 >= features.size() || features[i+1].featureName == "5UTR")
	    features[i].type = FIVEPRIMEUTRINITIAL;
	  else
	    features[i].type = FIVEPRIMEUTRSINGLE;
	}
	else if (i+1 >= features.size() || features[i+1].featureName == "5UTR")
	  features[i].type = FIVEPRIMEUTRINTERNAL;
	else
	  features[i].type = FIVEPRIMEUTRTERMINAL;
      }
      else {
	if (i == features.size() - 1) {
	  if (i==0 || features[i-1].featureName == "5UTR")
	    features[i].type = FIVEPRIMEUTRINITIAL;
	  else
	    features[i].type = FIVEPRIMEUTRSINGLE;
	}
	else {
	  if (i==0 || features[i-1].featureName == "5UTR")
	    features[i].type = FIVEPRIMEUTRINTERNAL;
	  else
	    features[i].type = FIVEPRIMEUTRTERMINAL;
	}
      }
    }
    else if (features[i].featureName == "3UTR") {
      features[i].type = THREEPRIMEUTR;
    }
    else if (features[i].featureName == "CDS") {
      if (strand == PLUS) {
	if (i > 0 && i < features.size() - 1 
	    && features[i-1].featureName == "start_codon" 
	    && features[i+1].featureName == "stop_codon") {
	  features[i].type = SINGLE;
	}
	else if (i > 0 && features[i-1].featureName == "start_codon") {
	  features[i].type = INITIAL;
	}
	else if (i < features.size() - 1 && features[i+1].featureName == "stop_codon") {
	  features[i].type = TERMINAL;
	}
	else {
	  features[i].type = INTERNAL;
	}
      }
      else {
	if (i > 0 && i < features.size() - 1 
	    && features[i-1].featureName == "stop_codon" 
	    && features[i+1].featureName == "start_codon") {
	  features[i].type = SINGLE;
	}
	else if (i > 0 && features[i-1].featureName == "stop_codon") {
	  features[i].type = TERMINAL;
	}
	else if (i < features.size() - 1 && features[i+1].featureName == "start_codon") {
	  features[i].type = INITIAL;
	}
	else {
	  features[i].type = INTERNAL;
	}
      }
    }
    else if (features[i].featureName == "start_codon" || 
	     features[i].featureName == "stop_codon") {
      //no exon types for these features
    }
    else {
      cerr << "No type for " << features[i].featureName << endl;
      fatalError ("Unrecognized feature name in GTF");
    }
  }
}


bool GTFTranscript::operator==(GTFTranscript& other) {
  if (id != other.id || start != other.start || end != other.end ||
      features.size() != other.features.size())
    return false;

  for (int i=0; i<features.size(); i++) {
    if (features[i] != other.features[i])
      return false;
  }

  return true;
}

bool GTFTranscript::operator!=(GTFTranscript& other) {
  return !operator==(other);
}

bool GTFTranscript::operator<(const GTFTranscript& other) const {
  return start < other.start;
}
