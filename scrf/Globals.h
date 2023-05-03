/* Global definitions.  Some data types, constants, and utility
   functions */

#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <iostream>
#include <math.h>
#include <cassert>
#include <new>
#include <map>

#ifdef QUAD_PRECISION
#include <qd/qd.h>
#endif

const int STRLEN = 1024;  /* the maximum length for short strings */

/* always use the std namespace */
using namespace std;

/* custom types that can be changed to balance between
   precision or max number of objects and memory usage */
typedef char stateid_t;
typedef int pos_t;
#ifdef QUAD_PRECISION
  typedef dd_real weight_t;
#else
  typedef double weight_t;
#endif
#define MPI_WEIGHT_T MPI_DOUBLE

typedef double small_weight_t;
const weight_t INF = 2e20;
const weight_t LOG_ZERO = -2e20;
const weight_t LOG_ZERO_REDUCED = -1e20;
const weight_t TINY = 1e-100;
enum sfsclass_t {MASKING, DNAKMER, DNAKMERPAIR, DNAKMERARRAY, DNAKMERPAIRARRAY, 
		 SVM, ESTPOSITION, ESTTRANSITION};
enum strand_t {PLUS, MINUS, NONE};
enum exon_t {INITIAL=0, INTERNAL, TERMINAL, SINGLE, 
	     FIVEPRIMEUTRINITIAL, FIVEPRIMEUTRINTERNAL, FIVEPRIMEUTRTERMINAL, FIVEPRIMEUTRSINGLE,
	     THREEPRIMEUTR};
static const string exonTypeNames[] = {"Initial", "Internal", "Terminal", "Single", 
				       "FivePrimeUTRInitial", "FivePrimeUTRInternal", "FivePrimeUTRTerminal", 
				       "FivePrimeUTRSingle", "ThreePrimeUTR"};
static const int numExonTypes = 9;
enum stat_t {CORRECT, PREDICTED, ANNOTATED};
static const int numStatTypes = 3;
enum kmer_requirement_t {NO_KMER_REQUIREMENT, TO_START, FROM_END};

/* when we hit an error, just quit  */
static void fatalError(string errorString) {
  cerr << "Fatal error: " << errorString << endl;
  exit(-1);
}

static inline bool** initializeNewBoolMatrix(int rows, int cols, bool value) {
  bool** matrix = new bool*[rows];
  bool** row_ptr_end = matrix + rows;
  for (bool** row_ptr=matrix; row_ptr != row_ptr_end; row_ptr++) {
    *row_ptr = new bool[cols];
    bool* col_ptr_end = *row_ptr + cols;
    for (bool* col_ptr=*row_ptr; col_ptr != col_ptr_end; col_ptr++)
      *col_ptr = value;
  }
  return matrix;
}

static inline weight_t** initializeNewWeightMatrix(int rows, int cols, weight_t value) {
  weight_t** matrix = new weight_t*[rows];
  weight_t** row_ptr;
  weight_t** row_ptr_end = matrix + rows;
  weight_t* col_ptr;
  weight_t* col_ptr_end;

  for (row_ptr=matrix; row_ptr != row_ptr_end; row_ptr++) {
    *row_ptr = new weight_t[cols];
    col_ptr_end = *row_ptr + cols;
    for (col_ptr=*row_ptr; col_ptr != col_ptr_end; col_ptr++)
      *col_ptr = value;
  }
  return matrix;
}

static inline void freeBoolMatrix(bool** matrix, int rows, int cols) {
  bool** row_ptr_end = matrix + rows;
  for (bool** row_ptr = matrix; row_ptr != row_ptr_end; row_ptr++)
    delete[] *row_ptr;
  delete[] matrix;
}

static inline void freeWeightMatrix(weight_t** matrix, int rows, int cols) {
  weight_t** row_ptr_end = matrix + rows;
  for (weight_t** row_ptr = matrix; row_ptr != row_ptr_end; row_ptr++)
    delete[] *row_ptr;
  delete[] matrix;
}

// rows and cols are the number of rows and cols in the input matrix
static inline void transposeBoolMatrix(bool** input, bool** output, int rows, int cols) {
  bool** out_row_ptr;
  int out_col = 0;

  bool** in_row_ptr_end = input + rows; 
  for (bool** in_row_ptr = input; in_row_ptr < in_row_ptr_end; in_row_ptr++) {
    bool* in_col_ptr_end = *in_row_ptr + cols;
    out_row_ptr = output;
    for (bool* in_col_ptr=*in_row_ptr; in_col_ptr != in_col_ptr_end; in_col_ptr++) {
      *(*out_row_ptr + out_col) = *in_col_ptr;
      out_row_ptr++;
    }
    out_col++;
  }
}

static inline void transposeWeightMatrix(weight_t** input, weight_t** output, int rows, int cols) {
  weight_t** out_row_ptr;
  int out_col = 0;

  weight_t** in_row_ptr_end = input + rows; 
  for (weight_t** in_row_ptr = input; in_row_ptr < in_row_ptr_end; in_row_ptr++) {
    weight_t* in_col_ptr_end = *in_row_ptr + cols;
    out_row_ptr = output;
    for (weight_t* in_col_ptr=*in_row_ptr; in_col_ptr != in_col_ptr_end; in_col_ptr++) {
      *(*out_row_ptr + out_col) = *in_col_ptr;
      out_row_ptr++;
    }
    out_col++;
  }
}

/* stuff for working with DNA and amino acid sequences */
const int DNA_CHARS = 7;

const int DNA_INDEX_COEFF[] = {1, 7, 49, 343, 2401, 16807, 117649, 823543};

const char INDEX_TO_DNA[] = {'A', 'C', 'G', 'T', '_', '.', 'N'};

const char DNA_TO_INDEX[] = {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 0, 6, 1, 6, 6, 6, 2, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4,
			     6, 0, 6, 1, 6, 6, 6, 2, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};

const char DNA_INDEX_COMPLEMENT[] = {3, 2, 1, 0, 4, 5, 6};

const int EST_CHARS = 5;
const int EST_INDEX_COEFF[] = {1, 5, 25, 125, 625};
const char INDEX_TO_EST[] = {'N', 'U', 'S', 'I', 'C'};
const char EST_TO_INDEX[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static void kmerIndexToString(int k, int kmerIndex, string& s) {
  s = "";
  for (int i=0; i<k; i++) {
    s += INDEX_TO_DNA[kmerIndex / DNA_INDEX_COEFF[k - i - 1]];
    kmerIndex -= DNA_INDEX_COEFF[k - i - 1] * (kmerIndex / DNA_INDEX_COEFF[k - i - 1]);
  }
  assert(kmerIndex == 0);
}

/* for any exon, given the 5' overhang and length % 3, returns
   the 3' overhang */
const int TOH_LOOKUP[3][3] = { {0, 1, 2}, {2, 0, 1}, {1, 2, 0} };

#endif
