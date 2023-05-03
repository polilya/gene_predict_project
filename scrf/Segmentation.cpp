#include "Segmentation.h"
#include <fstream>

Segmentation::Segmentation() {
}

void Segmentation::addSegment(stateid_t state, pos_t start, pos_t end) {
  if (state == -1) fatalError ("Unknown state in addSegment");

  Segment segment;
  segment.state = state;
  segment.start = start - 1; // argument is GTF coordinates, but we want 0-based coordinates
  segment.end = end - 1; //argument is GTF coordinates, but we wnat 0-based coordinates
  segments.push_back(segment);
}

Segmentation::Segmentation(GTF& gtf, SemiCRF& scrf) {
  length = scrf.alignment->length;

  gtf.discardOverlaps();     //overlapping segments not allowed
  gtf.discardShortExons(3);  //minimum exon length 3

  //create a vector of all transcripts
  vector<GTFTranscript> transcripts;
  for (map<string, GTFTranscript>::iterator iter = gtf.transcripts.begin();
       iter != gtf.transcripts.end(); ++iter) {
    transcripts.push_back(iter->second);
  }

  //sort transcripts by start coordinate
  sort(transcripts.begin(), transcripts.end());

  if (gtf.transcripts.size() == 0) {
    /* nothing annotated, assume all intergenic */
    addSegment(scrf.getStateID("Intergenic"), 1, scrf.alignment->length);
  }
  else {
    /* go through all transcripts and create appropriate segments */
    pos_t nextPosition = 1;
    for (int i=0; i<transcripts.size(); i++) {
      int codonPosition;  //keeps track of what codon position we're on for the regular CRF case
      if (transcripts[i].strand == PLUS)
	codonPosition = 0;
      else if (transcripts[i].strand == MINUS)
	codonPosition = 2;

      //add intergenic region before this transcript
      pos_t intergenicEnd = transcripts[i].start - 1;
      if (transcripts[i].strand == MINUS)
	intergenicEnd += 3;  //stop codon at start of transcript labeled intergenic
      addSegment(scrf.getStateID("Intergenic"), nextPosition, intergenicEnd);
      for (int j=0; j<transcripts[i].features.size(); j++) {
	if (transcripts[i].features[j].featureName == "5UTR") {
	  string stateName = exonTypeNames[transcripts[i].features[j].type];
	  if (transcripts[i].strand == MINUS)
	    stateName = stateName + "-";
	  addSegment(scrf.getStateID(stateName), transcripts[i].features[j].start,
		     transcripts[i].features[j].end);
	}
	else if (transcripts[i].features[j].featureName == "3UTR") {
	  if (transcripts[i].strand == PLUS)
	    addSegment(scrf.getStateID("ThreePrimeUTR"), transcripts[i].features[j].start,
		       transcripts[i].features[j].end);
	  else if (transcripts[i].strand == MINUS)
	    addSegment(scrf.getStateID("ThreePrimeUTR-"), transcripts[i].features[j].start,
		       transcripts[i].features[j].end);
	}
	else if (transcripts[i].features[j].featureName == "CDS") {
	  for (pos_t k=transcripts[i].features[j].start; k <= transcripts[i].features[j].end; k++) {
	    string stateName;
	    if (transcripts[i].features[j].type == INITIAL)
	      stateName = "InitialCDS";
	    else if (transcripts[i].features[j].type == INTERNAL)
	      stateName = "InternalCDS";
	    else if (transcripts[i].features[j].type == TERMINAL)
	      stateName = "TerminalCDS";
	    else if (transcripts[i].features[j].type == SINGLE)
	      stateName = "SingleCDS";
	    stateName += (char)codonPosition + 48;
	    if (transcripts[i].strand == PLUS) {
	      codonPosition++;
	      if (codonPosition > 2)
		codonPosition = 0;
	    }
	    else if (transcripts[i].strand == MINUS) {
	      codonPosition--;
	      if (codonPosition < 0)
		codonPosition = 2;
	      stateName += "-";
	    }

	    //special cases for the two most 5' bases of a start codon
	    if ((stateName == "InitialCDS0"||stateName == "SingleCDS0") && k == transcripts[i].features[j].start) {
	      addSegment(scrf.getStateID("StartCodonCDS0"), k, k);
	    }
	    else if ((stateName == "InitialCDS1"||stateName == "SingleCDS1") && k == transcripts[i].features[j].start + 1) {
	      addSegment(scrf.getStateID("StartCodonCDS1"), k, k);
	    }
	    else if ((stateName == "InitialCDS1-"||stateName == "SingleCDS1-") && k == transcripts[i].features[j].end - 1) {
	      addSegment(scrf.getStateID("StartCodonCDS1-"), k, k);
	    }
	    else if ((stateName == "InitialCDS0-"||stateName == "SingleCDS0-") && k == transcripts[i].features[j].end) {
	      addSegment(scrf.getStateID("StartCodonCDS0-"), k, k);
	    }
	    else {
	      addSegment(scrf.getStateID(stateName), k, k);
	    }
	  }
	}
	else {
	  continue;
	}

	nextPosition = transcripts[i].features[j].end + 1;
	if (j != transcripts[i].features.size() - 1 &&
	    nextPosition < transcripts[i].features[j+1].start) {
	  /* add an intron */
	  string intronName;
	  if (transcripts[i].features[j].featureName == "5UTR") {
	    if (transcripts[i].strand == PLUS)
	      intronName = "FivePrimeUTRIntron";
	    else
	      intronName = "FivePrimeUTRIntron-";
	  }
	  else if (scrf.states[segments.back().state].name == "ThreePrimeUTR")
	    intronName = "ThreePrimeUTRIntron";
	  else if (scrf.states[segments.back().state].name == "ThreePrimeUTR-")
	    intronName = "ThreePrimeUTRIntron-";
	  else {
	    /* coding sequence intron */
	    int phase;
	    if (scrf.numberOfImplicitStates == scrf.states.size())
	      phase = codonPosition;
	    else
	      phase = scrf.states[segments.back().state].phase;
	    intronName = "Intron";
	    intronName += (char)phase + 48;
	    if (transcripts[i].strand == MINUS)
	      intronName += "-";
	  }
	  addSegment(scrf.getStateID(intronName), nextPosition, 
		     transcripts[i].features[j+1].start - 1);
	  nextPosition = transcripts[i].features[j+1].start;
	}
      }
    }
    //add final intergenic region
    if (nextPosition <= scrf.alignment->length - 1)
      addSegment(scrf.getStateID("Intergenic"), nextPosition, scrf.alignment->length);
  }

  setLabels();
  //print(cerr, scrf);
}

void Segmentation::write(string outFilename, SemiCRF& scrf) {
  ofstream outStream(outFilename.c_str());
  print(outStream, scrf);
}

void Segmentation::print(ostream& outStream, SemiCRF& scrf) {

  for (int i=0; i<segments.size(); i++) {
    outStream << scrf.states[segments[i].state].name << "\t" << segments[i].start <<
      "\t" << segments[i].end << endl;
  }
}

bool Segmentation::setLabels() {
  labels.resize(length);

  for (pos_t i=0; i<length; i++)
    labels[i] = -1;  /* initialize to "unknown" */
  
  for (int i=0; i<segments.size(); i++) {
    for (pos_t j=segments[i].start; j<=segments[i].end; j++)
      labels[j] = segments[i].state;
  }
  
  bool success = true;

  for (pos_t i=0; i<length; i++) {
    if (labels[i] == -1) {
      cerr << "Label " << i << " was not set on Segmentation..." << endl;
      success = false;
    }
  }

  return success;
}


bool Segmentation::checkLegality(SemiCRF& scrf) {
  /* assume precomputeForDP has been called to set allowedTransitions */
  bool legal = true;

  setLabels();
  AlignmentSequence* alignment = scrf.alignment;

  if (length != alignment->length) {
    cerr << "Segmentation length did not match alignment length" << endl;
    legal = false;
  }

  for (pos_t j=1; j<alignment->length; j++) {
    stateid_t fromState = labels[j-1];
    stateid_t toState = labels[j];

    bool found = false;
    bool allowed = false;
    for (int t=0; t<scrf.transitions.size(); t++) {
      if (scrf.transitions[t].fromState == fromState && scrf.transitions[t].toState == toState) {
	found = true;
	if (scrf.dpMatrix.allowedTransition[j][t])
	  allowed = true;
      }
    }

    if (!found)
      cerr << "No transition found for labels " << scrf.states[fromState].name << ", " 
	   << scrf.states[toState].name << " at positions " << (j-1) << ", " << j << endl;
    if (found && !allowed)
      cerr << "Disallowed transition for labels " << scrf.states[fromState].name << ", " 
	   << scrf.states[toState].name << " at positions " << (j-1) << ", " << j << endl;

    if (!allowed)
      legal = false;
  }

  return legal;
}
