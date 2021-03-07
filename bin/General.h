

#ifndef _Off_Spotter_General_h
#define _Off_Spotter_General_h
//General shared memory variables, they will be set depending on the species
unsigned int SIZESKEY,HASHKEY,DATAKEY;

//General limits
#define MAX_CHROM_NUM 25        /* Maximum number of chromosomes (can be larger than actual number, eidt only if using larger genome)*/
#define KEY_NUM 4294967296      /* Number of ALL POSSIBLE different 16mers. This number is fixed for all genomes */
#define MAX_CHROM_LEN 268435455     /* Maximum chromosome length. Chromosome 1 of the human genome hg19 has length 249250622. This number uses 28 bits thus this max chrom len. Refer to POSBITS too. */
#define MAXPAMLENGTH 28             /* Maximum length of gRNA + PAM, edit when adding PAMs if new PAM is longer */
#define MYTYPE unsigned long        /* Renaming unsigned long, it is possible that we'll need something bigger in the future */

//Number of bits assigned per information type
#define LASTNTBITS 2        /* Each nucleotide from the last 4 nucleotides of a 20mer is represented by 2 bits 00=A,01=C,10=G,11=T */
#define EXTENSIONBITS 8     /* There are 4 nucleotides extra in a 20mer than a 16mer, and each counts for 2 bits */
#define STRANDBITS 1        /* Strand 0=+, 1=- */
#define CHROMBITS 5         /* Chromosome from 00000=0 to 11000=24 */
#define POSBITS 28          /* Position of hit within the chromosome. This is designed to fit a bit more that chr1 of hg19. If this is changed MAX_CHROM_LEN needs to be adjusted. */

/******************************************************* DIAGRAM OF EACH TABLE B ENTRY *******************************************/
/* ----------------------------------------------------------------------------------------------------------------------------- */
/* | Currently Unused | Last 3 Ns (or Rs) of PAM | PAM ID | Strand | First N of PAM | Chrom | Position in chrom | Last 4 Bases | */
/* ----------------------------------------------------------------------------------------------------------------------------- */
/* 64                 52                         46       44       43               41      36                  8              1 */
/*********************************************************************************************************************************/
//Shifting bits and masks. Used to recover the values above from a large number (that is from a table B entry)
//Each shift number value essentially is the number of bits preceeding a specific type of information (eg. chromosome = 36) and the
//mask is 0s 1s 0s and the 1s are placed only where the bits corresponding to a type of information is (eg chromosome = 23 0s, 5 1s, 36 0s)
//To recover the values we first apply one of the masks and then shift the number to the begiining of the table B entry.
//SHIFTS
#define POSSHIFTBITS 8
#define LASTNTSHIFTBITS 41
#define LAST3NTSHIFTBITS 46
#define PAMSHIFTBITS 44
#define CHROMSHIFTBITS 36
//MASKS
const long mask = 0xFFFFFFFFFFFFFFFF;       /* General purpose for TABLE B everything is 1. 64 1s */
const long keymask = 0xFFFFFFFF;            /* General purpose for TABLE A everything is 1. 32 1s */
const long extensionMask = 255;             /* Shift to get the last 4 bases. 56 0s, 8 1s */
const long posMask = 0xFFFFFFF00;           /* Shift mask to get the position. 28 0s, 28 1s, 8 0s */
const long chromMask = 0x1F000000000;       /* Shift to get the chromosome number. 23 0s, 5 1s, 36 0s */
const long lastntMask = 0x60000000000;      /* Shift to get the first N of the PAM. 21 0s, 2 1s, 41 0s */
const long strandMask = 0x80000000000;      /* Shift to get strand. 20 0s, 1 1s, 43 0s */
const long PAMmask = 0x300000000000;        /* Shift to get PAM. 18 0s, 2 1s, 44 0s */
const long last3ntMask = 0xC00000000000;    /* Shift to get the first from the right N of the PAM - The others are gotten recursively by shifting by 2 and then applying the mask. 16 0s, 2 1s, 46 0s */

//CONVERSION TABLES
//Every time a nucleotide needs to be converted eg. from C to G, we use the following tables.
//This can be done by finding it's ASCII value of the character and use these conversion tables to get a different character.
//Accordingly to get the char out of the int the reverse scheme is used, A=0,C=1,G=2,T=3
#define ASCIIA ('A')        /* The point of reference of all conversions for ASCII characters, the character A */
int cTOint [26] = {0,5,1,5,5,5,2,5,5,5,5,5,5,5,5,5,5,5,5,3,5,5,5,5,5,5};        /* This table assigns the value 0 to A, 1 to C, 2 to G and 3 to T to get the formula described above. 5 is the error number */
int cTOintREV [26] = {3,5,2,5,5,5,1,5,5,5,5,5,5,5,5,5,5,5,5,0,5,5,5,5,5,5};     /* This table assigns the value of the complement of each base, 3 to A (=T) etc */
char cTOcREV [26] = {'T','X','G','X','X','X','C','X','X','X','X','X','X','X','X','X','X','X','X','A','X','X','X','X','X','X'};; /* This table returns the reverse complement of a base */
char cTOnextC [26] = {'C','X','G','X','X','X','T','X','X','X','X','X','X','X','X','X','X','X','X','A','X','X','X','X','X','X'}; /* This table is used to circularly get the next nucleotide according to the rule A -> C -> G -> T -> A */
char cTOsmallC [26] = {'a','X','c','X','X','X','g','X','X','X','X','X','X','X','X','X','X','X','X','t','X','X','X','X','X','X'}; /* This table is used to get the corresponding small letter out of the capital one */
char intTOc [5] = {'A','C','G','T','X'};    /* This table returns the nucleotide that corresponds to the int */
char intTOpamID [4] = {'G', 'C', 'R', 'A'}; /* This table transforms an integer PAM ID to a char one. 0 = G = NGG, 1 = C = NNNNACA, 2 = R = NNGRRT, 3 = A = NAG */

#endif
