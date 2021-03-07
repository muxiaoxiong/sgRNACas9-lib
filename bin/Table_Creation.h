

#ifndef _Off_Spotter_filecreation_h
#define _Off_Spotter_filecreation_h

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "General.h"

using namespace std;

//FUNCTIONS
//To create binary files that hold table A and table B we need to first count all hits, then allocate the space, populate the tables and then
//write each on a separate file. Table A is the index file and table B the data.
//The files are created for convinience. If you are working with the same genome as the one used in the file creation last time you don't
//have to re-use this function, go to loading directly.
//If you want to add more or some of the PAMs, you need to create more such functions and then edit the main of Off-Spotter_Table_Creation.cpp
//as well as the Off-Spotter_Results files. Keep in mind that the total number of hits of all PAMs together should fit in Table B.
//Those functions count the number of hits on the input genome of a specific PAM for ALL possible gRNAs.
unsigned int countNNG(const char *filename, char PAM, unsigned int *currentPos);
int countNNNNACA(const char *filename, unsigned int *currentPos);
int countNNGRRT(const char *filename, unsigned int *currentPos);
//Those functions store the details of each hit on the input genome of a specific PAM and ALL possible gRNAs
void resultsNNG(const char *filename, char PAM, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID);
void resultsNNNNACA(const char *filename, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID);
void resultsNNGRRT(const char *filename, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID);

#endif
