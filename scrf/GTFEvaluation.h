#ifndef GTFEVALUATION_H
#define GTFEVALUATION_H

#include "Globals.h"
#include "GTF.h"
#include <iostream>

class GTFEvaluation {
 public:
  vector<int> transcriptStats;
  vector<int> combinedExonStats;
  vector<int> cdsPositionalStats;
  vector<vector<int> > individualExonStats;
  vector<vector<int> > exonTypePositionalStats;
  
  GTFEvaluation();
  GTFEvaluation(GTF& annotation, GTF& prediction);
  void pack(vector<int>& packed);
  void unpack(vector<int>& packed);
  void operator+=(GTFEvaluation& other);
  friend ostream& operator<<(ostream& os, GTFEvaluation& gtfEval);

 private:
  void init();  
};

#endif
