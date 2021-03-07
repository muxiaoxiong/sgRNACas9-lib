
#include <sys/ipc.h>
#include <sys/shm.h>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <unistd.h>

#include "Shared_Memory_IDs.h"
#include "General.h"

using namespace std;


int main (int argc, char *argv[]){
    char *usage = (char *) "\n\nUsage: %s -g <genome_name>\n\tPlease type hg38 or GRCh38 for human hg38\n\tPlease type hg19 or GRCh37 for human hg19\n\tPlease type mm10 or GRCm38 for mouse mm10\n\tPlease type w303 for yeast strain w303\n\n";
    //Parse input
    int a;
    char *genome = NULL;
    while ((a = getopt (argc, argv, "g:")) != -1){
        if(a == 'g'){
            genome=optarg;
            //Set the memory keys depending on the genome
            if (strcmp(genome,"hg19")==0 || strcmp(genome,"GRCh37")==0){
                HASHKEY=HASHKEYHUMAN;
                DATAKEY=DATAKEYHUMAN;
                SIZESKEY=SIZESKEYHUMAN;
            }
            else if (strcmp(genome,"mm10")==0 || strcmp(genome,"GRCm38")==0){
                HASHKEY=HASHKEYMOUSE;
                DATAKEY=DATAKEYMOUSE;
                SIZESKEY=SIZESKEYMOUSE;
            }
            else if (strcmp(genome,"hg38")==0 || strcmp(genome,"GRCh38")==0){
                HASHKEY=HASHKEYHUMAN2;
                DATAKEY=DATAKEYHUMAN2;
                SIZESKEY=SIZESKEYHUMAN2;
            }
            else if (strcmp(genome,"w303")==0){
                HASHKEY=HASHKEYYEAST;
                DATAKEY=DATAKEYYEAST;
                SIZESKEY=SIZESKEYYEAST;
            }
            else{
                fprintf(stderr, "\nUknown <genome_name> %s\n", genome);
                fprintf(stderr, usage, argv[0]);
                fflush(stderr);
                return(-1);
            }
        }
        else if (a == '?'){
            if (optopt == 'c')
                fprintf (stderr, "\nOption -%c requires an argument\n", optopt);
            else if (isprint (optopt))
                fprintf (stderr, "\nUnknown option `-%c'\n", optopt);
            else
                fprintf (stderr, "\nUnknown option character `\\x%x'\n", optopt);
            fprintf(stderr, usage, argv[0]);
            fflush(stderr);
            return -1;
        }
        else{
            fprintf(stderr, usage, argv[0]);
            fflush(stderr);
            return(1);
        }
    }
    if (genome == NULL){
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }

    //Get the size of table B
    int sizes_ID;
    if ((sizes_ID = shmget(SIZESKEY,sizeof(unsigned int), 0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s are loaded on memory\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    unsigned int *sizes;
    if ((sizes = (unsigned int *) shmat(sizes_ID,NULL,0))==  (void *)-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s are loaded on memory\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }

    //Dettach table A
    cout << "\nDetaching Table A" << endl;
    int ID;
    if (( ID = shmget(HASHKEY,KEY_NUM*sizeof(unsigned int), 0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s are loaded on memory\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    if (shmctl(ID, IPC_RMID, 0) == -1){
        fprintf(stderr, "\nAn error occured while deleting memory segment of %s\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }

    //Dettach table B
    cout << "\nDetaching Table B" << endl;
    if ((ID = shmget(DATAKEY,sizes[0]*sizeof(MYTYPE),0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s are loaded on memory\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    if (shmctl(ID, IPC_RMID, 0) == -1){
        fprintf(stderr, "\nAn error occured while deleting memory segment of %s\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }

    //Dettach the size variable
    if (shmctl(sizes_ID, IPC_RMID, 0)){
        fprintf(stderr, "\nAn error occured while deleting memory segment of %s\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }

    cout << "\nDone!" << endl;

    return 1;
}
