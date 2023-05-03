/* A sementation, or labeling, of the sequence */

#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include "GTF.h"
#include "Segment.h"
#include "Globals.h"
#include "State.h"
#include "SemiCRF.h"
#include <string>
#include <vector>
#include <algorithm>

class GTF;  /* forward declarations */
class SemiCRF;

class Segmentation {
 public:
  vector<Segment> segments;
  vector<stateid_t> labels; /* state ID at each position.  -1 means "unknown" */
  pos_t length;

  Segmentation();
  Segmentation(GTF& gtf, SemiCRF& scrf);
  void addSegment(stateid_t state, pos_t start, pos_t end);
  void write(string outFilename, SemiCRF& scrf);
  void print(ostream& outStream, SemiCRF& scrf);
  bool setLabels(); /* set the labels based on the segments */
  bool checkLegality(SemiCRF& scrf);  //are all transitions allowed based on the sequence?
};

#endif
