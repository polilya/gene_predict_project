#ifndef GTF_H
#define GTF_H

#include "SemiCRF.h"
#include "Segmentation.h"
#include "GTFTranscript.h"
#include <vector>
#include <string>
#include <map>
#include <iostream>

class Segmentation;  /* forward declaration */
class SemiCRF;       /* forward declaration */

class GTF {
 public:
  map<string, GTFTranscript> transcripts;

  GTF(string gtfFilename);
  GTF(Segmentation& segmentation, SemiCRF& scrf);
  void write(string outFilename);
  void print(ostream& outStream);
  void discardOverlaps(); /* randomly discards overlapping transcripts so
			     that final GTF contains no overlaps */
  void discardShortExons(int minLength); /* discards any transcripts containing 
					    exons shorter than minLength */
  void removeFeatureType(string featureName);
  
 private:
  /* helper functions */
  void addStartCodon(pos_t start, strand_t strand, string transcriptID);
  void addExon(pos_t start, pos_t end, weight_t score, strand_t strand, string transcriptID);
  void addStopCodon(pos_t end, strand_t strand, string transcriptID);
  void addFivePrimeUTR(pos_t start, pos_t end, strand_t strand, string transcriptID);
  void addThreePrimeUTR(pos_t start, pos_t end, strand_t strand, string transcriptID);
};

#endif
