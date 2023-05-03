#ifndef GTFFEATURE_H
#define GTFFEATURE_H

#include "Globals.h"
#include <string>
#include <fstream>

class GTFFeature {
 public:
  string seqName;
  string source;
  string featureName;
  pos_t start;
  pos_t end;
  weight_t score;
  strand_t strand;
  int foh;  /* five-prime overhang, i.e., how many bases
	       on the five-prime end before the first complete
	       codon */
  string geneID;
  string transcriptID;
  bool scoreDefined;
  exon_t type;

  void read(ifstream& gtfStream);
  void print(ostream& outStream);

  bool operator==(const GTFFeature& other);
  bool operator!=(const GTFFeature& other);
  bool operator<(const GTFFeature& other) const;
};

#endif
