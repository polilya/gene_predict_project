#include "GTF.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <time.h>

GTF::GTF(string gtfFilename) {
  ifstream gtfStream(gtfFilename.c_str());

  if (gtfStream.fail()) {
    cerr << "Could not open " << gtfFilename << endl;
    fatalError("GTF file open error");
  }

  while (! gtfStream.eof()) {
    GTFFeature feature;
    feature.read(gtfStream);
    if (! gtfStream.eof())
      transcripts[feature.transcriptID].features.push_back(feature);
  }

  /* sort the features in each transcript */
  for (map<string, GTFTranscript>::iterator iter = transcripts.begin();
       iter != transcripts.end(); ++iter) {
    sort(iter->second.features.begin(), iter->second.features.end());
  }

  /* set transcript IDs, coordinates, and strand */
   for (map<string, GTFTranscript>::iterator iter = transcripts.begin();
       iter != transcripts.end(); ++iter) {
     iter->second.setAttributesFromFeatures();
   }

  /* throw out any incomplete transcripts */
  /* eventually we will want to allow these to stay, but right now
     the code doesn't handle them */
  for (map<string, GTFTranscript>::iterator iter = transcripts.begin();
       iter != transcripts.end(); ) {
    GTFTranscript& transcript = iter->second;
    int numStartCodons = 0;
    int numStopCodons = 0;
    for (int i=0; i<transcript.features.size(); i++) {
      if (transcript.features[i].featureName == "start_codon")
	numStartCodons++;
      else if (transcript.features[i].featureName == "stop_codon")
	numStopCodons++;
    }
    if (numStartCodons != 1 || numStopCodons != 1) {
      //cerr << "Discarding incomplete transcript " << transcript.id << endl;
      transcripts.erase(iter++);
    }
    else
      ++iter;
  }

  /* set correct five-prime overhangs based on length in case they are missing or
     wrong */
  /* also set exon type information */
  for (map<string, GTFTranscript>::iterator iter = transcripts.begin();
       iter != transcripts.end(); ++iter) {
    iter->second.setFohs();
    iter->second.setExonTypes();
  }
}

GTF::GTF(Segmentation& segmentation, SemiCRF& scrf) {
  int transcriptNum = 0;
  string transcriptID = "CONTRAST.0";

  for (int i=0; i<segmentation.segments.size(); i++) {
    Segment& seg = segmentation.segments[i];
    string stateName = scrf.states[seg.state].name;
    if (stateName == "Intergenic") {
      transcriptNum++;
      ostringstream transcriptIDStream;
      transcriptIDStream << "CONTRAST." << transcriptNum;
      transcriptID = transcriptIDStream.str();
    }

    if (stateName == "FivePrimeUTR" || stateName == "FivePrimeUTR-" ||
	stateName == "FivePrimeUTRInitial" || stateName == "FivePrimeUTRInitial-" ||
	stateName == "FivePrimeUTRInternal" || stateName == "FivePrimeUTRInternal-" ||
	stateName == "FivePrimeUTRTerminal" || stateName == "FivePrimeUTRTerminal-" ||
	stateName == "FivePrimeUTRSingle" || stateName == "FivePrimeUTRSingle-" )
      addFivePrimeUTR(seg.start, seg.end, scrf.states[seg.state].strand, transcriptID);

    if (stateName == "ThreePrimeUTR" || stateName == "ThreePrimeUTR-")
      addThreePrimeUTR(seg.start, seg.end, scrf.states[seg.state].strand, transcriptID);

    if (stateName.find("CDS", 0) != string::npos) { 
      if (i > 0) {
	string prevName = scrf.states[segmentation.segments[i-1].state].name;
	if (prevName.substr(0,6) != "Intron") {
	  /* beginning new gene */
	  if (scrf.states[seg.state].strand == PLUS)
	    addStartCodon(seg.start, scrf.states[seg.state].strand, transcriptID);
	  else if (scrf.states[seg.state].strand == MINUS)
	    addStopCodon(seg.start-3, scrf.states[seg.state].strand, transcriptID);
	}
      }      
      
      pos_t cdsStart = seg.start;
      while (i < segmentation.segments.size() && scrf.states[segmentation.segments[i].state].name.find("CDS", 0) != string::npos)
	i++;
      i--;
      pos_t cdsEnd = segmentation.segments[i].end;

      weight_t score = -1;
      if (scrf.dpMatrix.alpha != NULL && scrf.dpMatrix.beta != NULL) {
	//compute exon posterior probability to use as GTF score
	
	weight_t exonWeight = 0;

	for (int j=cdsStart; j<=cdsEnd+1 && j<scrf.alignment->length; j++) {
	  stateid_t thisState = segmentation.labels[j];
	  
	  exonWeight += scrf.dpMatrix.stateSFSSums[j][thisState];
	  
	  if (j == 0)
	    exonWeight += scrf.states[thisState].startWeight;
	  else {
	    stateid_t prevState = segmentation.labels[j-1];
	    for (int k=0; k<scrf.states[thisState].transitionsTo.size(); k++) {
	      int transitionID = scrf.states[thisState].transitionsTo[k];
	      TransitionFeature& transition = scrf.transitions[transitionID];
	      if (transition.fromState == prevState && scrf.dpMatrix.allowedTransition[j][transitionID]) {
		exonWeight += transition.weight;
		
		weight_t transitionSFSSum = transition.typeID == -1 ? (weight_t)0 : 
		  scrf.dpMatrix.transitionSFSSums[j][2 * transition.typeID + transition.strand];
		exonWeight += transitionSFSSum;
	      }
	    }
	  }
	}
	
	if (cdsStart == 0 && cdsEnd == scrf.alignment->length-1) {
	  score = exp(exonWeight - scrf.computePartitionFromAlpha());
	}
	else if (cdsStart == 0) {
	  score = exp(exonWeight +
		      scrf.dpMatrix.beta[cdsEnd+1][segmentation.labels[cdsEnd+1]] - 
		      scrf.computePartitionFromAlpha());
	}
	else if (cdsEnd == scrf.alignment->length - 1) {
	  score = exp(scrf.dpMatrix.alpha[cdsStart-1][segmentation.labels[cdsStart-1]] +
		      exonWeight - 
		      scrf.computePartitionFromAlpha());
	}
	else {
	  score = exp(scrf.dpMatrix.alpha[cdsStart-1][segmentation.labels[cdsStart-1]] +
		      exonWeight +
		      scrf.dpMatrix.beta[cdsEnd+1][segmentation.labels[cdsEnd+1]] - 
		      scrf.computePartitionFromAlpha());
	}
      }
      
      addExon(cdsStart, cdsEnd, score, scrf.states[seg.state].strand, transcriptID);
      
      seg = segmentation.segments[i];
      stateName = scrf.states[seg.state].name;
      
      if (i+1 < segmentation.segments.size()) {
	string nextName = scrf.states[segmentation.segments[i+1].state].name;
	if (nextName.substr(0,6) != "Intron") {
	  /* ending gene */
	  if (scrf.states[seg.state].strand == PLUS)
	    addStopCodon(seg.end + 1, scrf.states[seg.state].strand, transcriptID);
	  else if (scrf.states[seg.state].strand == MINUS)
	    addStartCodon(seg.start - 2, scrf.states[seg.state].strand, transcriptID);
	}
      }
    }
  }

  /* sort the features in each transcript */
  for (map<string, GTFTranscript>::iterator iter = transcripts.begin();
       iter != transcripts.end(); ++iter) {
    sort(iter->second.features.begin(), iter->second.features.end());
  }

  /* set transcript IDs, coordinates, and strand */
  for (map<string, GTFTranscript>::iterator iter = transcripts.begin();
       iter != transcripts.end(); ++iter) {
    iter->second.setAttributesFromFeatures();
  }
  
  /* set correct five-prime overhangs based on length in case they are missing or
     wrong */
  /* also set exon type information */
  for (map<string, GTFTranscript>::iterator iter = transcripts.begin();
       iter != transcripts.end(); ++iter) {
    iter->second.setFohs();
    iter->second.setExonTypes();
  }
}
   
void GTF::addStartCodon(pos_t start, strand_t strand, string transcriptID) {
  GTFFeature startCodon;

  startCodon.seqName = "SomeSequence";
  startCodon.source = "CONTRAST";
  startCodon.featureName = "start_codon";
  startCodon.start = start + 1;   /* GTF starts at coordinate 1 */
  startCodon.end = start + 2 + 1;
  startCodon.score = (weight_t)0;
  startCodon.scoreDefined = false;
  startCodon.strand = strand;
  startCodon.foh = 0;
  startCodon.geneID = transcriptID;
  startCodon.transcriptID = transcriptID;

  transcripts[transcriptID].features.push_back(startCodon);
}

void GTF::addExon(pos_t start, pos_t end, weight_t score, strand_t strand, string transcriptID) {
  GTFFeature exon;

  exon.seqName = "SomeSequence";
  exon.source = "CONTRAST";
  exon.featureName = "CDS";
  exon.start = start + 1;   /* GTF starts at coordinate 1 */
  exon.end = end + 1;

  if (score < 0) {
    exon.score = (weight_t)0;
    exon.scoreDefined = false;
  }
  else {
    exon.score = score;
    exon.scoreDefined = true;
  }

  exon.strand = strand;
  exon.geneID = transcriptID;
  exon.transcriptID = transcriptID;

  transcripts[transcriptID].features.push_back(exon);
}

void GTF::addStopCodon(pos_t start, strand_t strand, string transcriptID) {
  GTFFeature stopCodon;

  stopCodon.seqName = "SomeSequence";
  stopCodon.source = "CONTRAST";
  stopCodon.featureName = "stop_codon";
  stopCodon.start = start + 1;   /* GTF starts at coordinate 1 */
  stopCodon.end = start + 2 + 1;
  stopCodon.score = (weight_t)0;
  stopCodon.scoreDefined = false;
  stopCodon.strand = strand;
  stopCodon.foh = 0;
  stopCodon.geneID = transcriptID;
  stopCodon.transcriptID = transcriptID;

  transcripts[transcriptID].features.push_back(stopCodon);
}

void GTF::addFivePrimeUTR(pos_t start, pos_t end, strand_t strand, string transcriptID) {
  GTFFeature utr;

  utr.seqName = "SomeSequence";
  utr.source = "CONTRAST";
  utr.featureName = "5UTR";
  utr.start = start + 1;   /* GTF starts at coordinate 1 */
  utr.end = end + 1;
  utr.score = (weight_t)0;
  utr.scoreDefined = false;
  utr.strand = strand;
  utr.foh = 0;
  utr.geneID = transcriptID;
  utr.transcriptID = transcriptID;

  transcripts[transcriptID].features.push_back(utr);
}

void GTF::addThreePrimeUTR(pos_t start, pos_t end, strand_t strand, string transcriptID) {
  GTFFeature utr;

  utr.seqName = "SomeSequence";
  utr.source = "CONTRAST";
  utr.featureName = "3UTR";
  utr.start = start + 1;   /* GTF starts at coordinate 1 */
  utr.end = end + 1;
  utr.score = (weight_t)0;
  utr.scoreDefined = false;
  utr.strand = strand;
  utr.foh = 0;
  utr.geneID = transcriptID;
  utr.transcriptID = transcriptID;

  transcripts[transcriptID].features.push_back(utr);
}

void GTF::write(string outFilename) {
  ofstream outStream(outFilename.c_str());
  print(outStream);
}

void GTF::print(ostream& outStream) {
  //create a vector of all transcripts
  vector<GTFTranscript> transcriptVector;
  for (map<string, GTFTranscript>::iterator iter = transcripts.begin();
       iter != transcripts.end(); ++iter) {
    transcriptVector.push_back(iter->second);
  }
  //sort transcripts by start coordinate
  sort(transcriptVector.begin(), transcriptVector.end());
  
  for (int i=0; i<transcriptVector.size(); i++)
    transcriptVector[i].print(outStream);
}

void GTF::removeFeatureType(string featureName) {
  for (map<string, GTFTranscript>::iterator i = transcripts.begin();
       i != transcripts.end(); i++) {
    vector<GTFFeature>& features = i->second.features;
    vector<GTFFeature>::iterator it = features.begin();
    while (it != features.end()) {
      if (it->featureName == featureName) {
	it = features.erase(it);
      }
      else {
	++it;
      }
    }
    i->second.setAttributesFromFeatures();  //start and stop may have changed
  }
}

void GTF::discardOverlaps() {
  /* assign all transcripts a unique random priority to determine which
     will be deleted in case of an overlap */
  map<string, int> transcriptPriority;
  for (map<string, GTFTranscript>::iterator i = transcripts.begin();
       i != transcripts.end(); i++) {
    bool duplicate;
    do {
      duplicate = false;
      transcriptPriority[i->first] = rand();
      for (map<string, int>::iterator j = transcriptPriority.begin();
	   j != transcriptPriority.end(); j++) {
	if (i->first != j->first && transcriptPriority[i->first] == j->second)
	  duplicate = true;
      }
    } while (duplicate);
  }

  /* Now check for overlaps and build up a list of transcripts to
     discard */
  set<string> discards;
  for (map<string, GTFTranscript>::iterator i = transcripts.begin();
       i != transcripts.end(); ++i) {
    for (map<string, GTFTranscript>::iterator j = transcripts.begin();
	 j != transcripts.end(); ++j) {
      if (i == j)
	continue;
      if (i->second.overlaps(j->second)) {
	/* we have two overlapping transcripts, mark the lower priority one
	 for deletion */
	if (transcriptPriority[i->first] < transcriptPriority[j->first])
	  discards.insert(i->first);
	else
	  discards.insert(j->first);
      }
    }
  }

  for (set<string>::iterator i = discards.begin(); i != discards.end(); i++) {
    //cerr << "Discarding overlapping transcript " << *i << endl;
    transcripts.erase(*i);
  }
}

void GTF::discardShortExons(int minLength) {
  set<string> discards;
  
  for (map<string, GTFTranscript>::iterator i = transcripts.begin();
       i != transcripts.end(); ++i) {
    bool discard = false;
    vector<GTFFeature>& features = i->second.features;
    for (int j=0; j<features.size(); j++) {
      int length = features[j].end - features[j].start + 1;
      if (features[j].featureName == "CDS" && length < 3)
	discard = true;
    }
    if (discard)
      discards.insert(i->first);
  }

  for (set<string>::iterator i = discards.begin(); i != discards.end(); i++) {
    //cerr << "Discarding transcript with short exon: " << *i << endl;
    transcripts.erase(*i);
  }
}
