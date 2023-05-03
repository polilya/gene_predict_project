#include <new>
#include <map>
#include <string>
#include <sstream>
#include "Globals.h"

#include "MemoryDebug.h"

#undef new

using namespace std;

bool mdTrackMemory = false;
map<void*, MemoryBlock> mdArrayBlocks;
map<void*, MemoryBlock> mdStandardBlocks;

string mdFile = "Unknown";
int mdLine = 0;

void showMemoryUsage() {
  mdTrackMemory = false;

  map<string, int> lineUsage;

  for (map<void*, MemoryBlock>::iterator iter = mdStandardBlocks.begin(); iter != mdStandardBlocks.end(); ++iter) {
    ostringstream lineStream;
    lineStream << iter->second.file << ":" << iter->second.line;
    lineUsage[lineStream.str()] += iter->second.size;
  }
  for (map<void*, MemoryBlock>::iterator iter = mdArrayBlocks.begin(); iter != mdArrayBlocks.end(); ++iter) {
    ostringstream lineStream;
    lineStream << iter->second.file << ":" << iter->second.line;
    lineUsage[lineStream.str()] += iter->second.size;
  }
  for (map<string, int>::iterator iter = lineUsage.begin(); iter != lineUsage.end(); ++iter) {
    cerr << iter->first << "\t" << (double)iter->second/1e6 << endl;
  }

  mdTrackMemory = true;
}

int mdSetLineAndFile(string file, int line) {
  cerr << "Setting file " << file << " and line " << line << endl;
  mdFile = file;
  mdLine = line;
  return 0;
}

void* mdNew(size_t size, map<void*, MemoryBlock>& mdMap) throw(bad_alloc) {

  void* ptr = malloc(size);

  if (mdTrackMemory) {
    mdTrackMemory = false;
    //if (size == 0) {
    //     cerr << "Allocation of size 0 called from " << mdFile << ":" << mdLine << endl;
    //}
    cerr << "Allocating block of size " << size << " from " << mdFile << endl;
    MemoryBlock block(mdFile, mdLine, size);
    mdMap[ptr] = block; 
    mdFile = "Unknown";
    mdLine = 0;
    mdTrackMemory = true;
  }

  return ptr;
}

void* mdDelete(void* ptr, map<void*, MemoryBlock>& mdMap) {

  if (mdTrackMemory) {
    mdTrackMemory = false;
    if (mdMap.find(ptr) == mdMap.end()) {
      cerr << "Double free or incorrect delete type" << endl;
      exit(0);
    }
    mdMap.erase(ptr);
    mdTrackMemory = true;
  }

  free(ptr);
}

void* operator new(size_t size) throw(bad_alloc) {
  return mdNew(size, mdStandardBlocks);
}

void* operator new[](size_t size) throw(bad_alloc) {
  return mdNew(size, mdArrayBlocks);
}

void operator delete(void* ptr) {
  mdDelete(ptr, mdStandardBlocks);
}

void operator delete[](void* ptr) {
  mdDelete(ptr, mdArrayBlocks);
} 

