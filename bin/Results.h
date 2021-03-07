
#ifndef __results_classes__
#define __results_classes__

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <cerrno>
#include <vector>
#include <string>
#include <sstream>

#include "General.h"
#include "Shared_Memory_IDs.h"

// This is an important number that can affect the annotation's correctness. Imagine that there are two entries in
// an annotation file e1 and e2. The coordinates of e1 are a and b and of e2 c and d, where a and c are the start
// coordinates. The annotation files are ordered by chromosome then by the start coordinate and then by the end
// coordinate. That means that a < c ALWAYS but it is possible that b > d. For example: e1: 1500-1800 and
// e2: 1510-1600. The result printing algorithm has to look through the annotation files for each result in order to
// find the correct annotations. To do that it loads the annotation file once and keeps an index table of the start
// and end line for each chromosome. Since the entries for each chromosome are ordered if b < d ALWAYS then we could
// do binary search to find the annotations to be printed. However now we need to do binary search on a range which
// is more complicated. Statistically we noticed that for the annotation files we used the number of entries between
// each e1 and e2 is small, say g. Thus, a simple way to solve the problem is to do binary search on the first
// coordinate and then move g entries towards the top. This results in O(n) time in the worst case where e1 is the
// first and e2 is the last entry which is the time that it would take if we didn't use binary search at all.
// However we think that it will speed up the process if used for any annotation file. If used it should always be
// equal to or larger of the biggest gap g. If not used it should be set to any number than the maximum number of
// entries for a chromosome in the annotation files, for example the length of the larger annotation file used.  If
// no annotation files are used it doesn't affect anything. Here we set it to max unsigned int to fit the general case.
#define MAX_BACKWARD_GAP 65535
#define MAX_BACKWARD_GAP_HG19 58
#define MAX_BACKWARD_GAP_MM10 106
#define MAX_BACKWARD_GAP_HG38 107

using namespace std;

//Table A
unsigned int *hashTable;
//Table B
MYTYPE *tableValues;
//Total length of table B
unsigned int total_len;

//Class for results. Each has the chromosome number (0-24), the strand (0 or 1), the in chromosome position and
//annotation info
class coords {
public:
   int chrom;
   int strand;
   unsigned int start;
   unsigned int end;
   std::string info;
};

//F U N C T I O N S
//Attaches to the already existing shared memory segments of table A, table B and size of table B.
void attach_shared_memory();
//Assigns the annotation info to the results found before they are returned
void attach_class_details(vector<string> *results, char* filepath, unsigned int backward_gap);
//Finds the results from shared memory that satisfy the query
int find_results (char * query,char *orig_query, int dots, int p1, int p2,int p3, int p4, int p5, char PAM, vector<string> *results);
//Runs the above function for all gRNAs that have at most "mismatches" mismatches from the query
void wildcard_find(char * query, char PAM,vector<string> *results, int mismatches);

#endif /* defined(_____0mer_results_classes__) */
