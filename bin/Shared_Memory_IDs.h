


// Header used by Load memory, Results and Dettach memory.

#ifndef _Off_Spotter_Shared_Memory_IDs_h
#define _Off_Spotter_Shared_Memory_IDs_h

// Those Ids are used by the shared memory functions to identify the shared memory segments
// across the programs.

/**********************************************HUMAN********************************************************/
// HASHKEY is the ID for Table A, the index table.
#define HASHKEYHUMAN 4231
// DATAKEY is the ID for Table B, the table that holds the results.
#define DATAKEYHUMAN 8675
//SIZESKEY holds the total number of entries of Table B.
#define SIZESKEYHUMAN 9000

/**********************************************MOUSE********************************************************/
// HASHKEY is the ID for Table A, the index table.
#define HASHKEYMOUSE 4232
// DATAKEY is the ID for Table B, the table that holds the results.
#define DATAKEYMOUSE 8676
//SIZESKEY holds the total number of entries of Table B.
#define SIZESKEYMOUSE 8999

/**********************************************HUMAN2*******************************************************/
// HASHKEY is the ID for Table A, the index table.
#define HASHKEYHUMAN2 4233
// DATAKEY is the ID for Table B, the table that holds the results.
#define DATAKEYHUMAN2 8677
//SIZESKEY holds the total number of entries of Table B.
#define SIZESKEYHUMAN2 8998

/**********************************************YEAST********************************************************/
// HASHKEY is the ID for Table A, the index table.
#define HASHKEYYEAST 4234
// DATAKEY is the ID for Table B, the table that holds the results.
#define DATAKEYYEAST 8678
//SIZESKEY holds the total number of entries of Table B.
#define SIZESKEYYEAST 8997

#endif
