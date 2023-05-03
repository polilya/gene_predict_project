#ifndef GTFTRANSCRIPT_H
#define GTFTRANSCRIPT_H

#include "Globals.h"
#include "GTFFeature.h"
#include <vector>
#include <string>
#include <fstream>

class GTFTranscript {
 public:
  string id;
  pos_t start;
  pos_t end;
  strand_t strand;
  vector<GTFFeature> features;

  void setAttributesFromFeatures(); //set ID, start, end, and strand based on features
  void print(ostream& outStream);
  bool overlaps(GTFTranscript& other);
  void setFohs();  /* set five-prime overhang information */
  void setExonTypes();

  bool operator==(GTFTranscript& other);
  bool operator!=(GTFTranscript& other);
  bool operator<(const GTFTranscript& other) const;
};

#endif
