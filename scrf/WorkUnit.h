#ifndef WORKUNIT_H
#define WORKUNIT_H

#include "Globals.h"
#include <time.h>

class WorkUnit {
 public:
  int id;
  string alignmentFile;
  string annotationFile;
  int assignedProcess;
  time_t startTime;
  time_t endTime;
  time_t elapsedTime;

  WorkUnit(int id, string alignmentFile, string annotationFile, int assignedProcess);
  bool operator<(const WorkUnit& other) const;
};

#endif
