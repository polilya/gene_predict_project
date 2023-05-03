#include "WorkUnit.h"

WorkUnit::WorkUnit(int id, string alignmentFile, string annotationFile, 
		   int assignedProcess) {
  this->id = id;
  this->alignmentFile = alignmentFile;
  this->annotationFile = annotationFile;
  this->assignedProcess = assignedProcess;
  elapsedTime = 0;
}

bool WorkUnit::operator<(const WorkUnit& other) const {
  return elapsedTime > other.elapsedTime;  /* sort from longest to shortest time */
}
