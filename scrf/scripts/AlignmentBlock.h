#include <vector>
#include <iostream>

using namespace std;

class AlignmentBlock {
 public:
  double score;
  unsigned long start;
  unsigned long length;
  vector<int> species;
  vector<char*> sequences;
  vector<unsigned long> compressedSequenceLengths;

  bool operator<(const AlignmentBlock& other) const;
  void addSequence(int species, string& sequence);
  void print(ostream& os);
  void compressBlock();
  void uncompressBlock();
};
