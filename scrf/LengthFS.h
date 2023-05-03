/* A set of features based on segment length */

#ifndef LENGTHFS_H
#define LENGTHFS_H

#include "FeatureSet.h"

class LengthFS : public FeatureSet {
 public:
  int maxLength;
  int binSize;    //group binSize lengths together as one parameter

  LengthFS(int maxLength, string parameterization);
  virtual void setValueToWeight();
};

#endif
