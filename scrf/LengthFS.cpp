#include "LengthFS.h"
#include <math.h>

LengthFS::LengthFS(int maxLength, string parameterization) {
  this->maxLength = maxLength;
  this->parameterization = parameterization;
  binSize = 10;

  int numValues = maxLength + 1;
  valueToParam.resize(numValues);
  valueToWeight.resize(numValues);
  valueToGlobalParamID.resize(numValues);

  if (parameterization == "Geometric") {
    for (int i=0; i<numValues; i++)
      valueToParam[i] = 0; /* not used */
  }
  else if (parameterization == "Standard") {
    /* bin size 10 */
    for (int i=1; i<numValues; i++)
      valueToParam[i] = (i-1) / binSize;
  }
  else
    fatalError("Unknown parameterization for LengthFS");
}

void LengthFS::setValueToWeight() {
  for (int i=0; i<valueToWeight.size(); i++) {
    if (parameterization == "Geometric") 
      valueToWeight[i] = paramToWeight[0] + (weight_t)(i-1) * paramToWeight[1];
    else
      valueToWeight[i] = paramToWeight[valueToParam[i]];
  }
}
