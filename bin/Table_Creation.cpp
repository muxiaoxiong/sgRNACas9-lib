
#include <getopt.h>
#include <unistd.h>

#include "Table_Creation.h"


using namespace std;


/***********************************************************************************************************/
/*******************************************   COUNT FUNCTIONS   *******************************************/
/***********************************************************************************************************/
/* Those functions count the amount of results per PAM per 16base key. The results of all the functions    */
/* kept in the same table (one counter per key), this table is currentPos. total_len is a counter that has */
/* the number of all results of all PAMs. Thus in the end of each of those functions we know how many hits */
/* there are in the genome per 16mer for this PAM and their total number. We need this info to allocate    */
/* enough space for table B and populate with the sum of hits of each individual 16mer for all PAMs table A*/
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/

/***********************************************************************************************************/
/* This function counts the appearances of ALL targets of ALL possible gRNAs followed by NGG or NAG on the */
/* filename argument. CurrentPos is keeping counts for each gRNA. It returns the total number.             */
/***********************************************************************************************************/
unsigned int countNNG(const char *filename, char PAM, unsigned int *currentPos){
    unsigned int i,j,end,size,key,total_num=0;
    char RevPAM;
    string chrom(MAX_CHROM_LEN,'\0');
    MYTYPE twobit,value,rest4;
    bool insert_entry;

    //Depending on the input PAM find reverse for 2nd N
    if (PAM=='G')
        RevPAM = 'C';
    else if (PAM == 'A')
        RevPAM = 'T';
    else
        return 0;

    //Open the input file for reading
    ifstream ifile (filename,std::ifstream::in);

    //Count appearances per 20mer
    //Forward
    while (!ifile.eof()){
        //header
        getline(ifile,chrom);
        if (ifile.eof()) break;
        //all lines up to next header
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 22){
            if(chrom[i+21] == PAM && chrom[i+22] == 'G'){
                insert_entry = true;
                //Make the key
                end = i + 16;
                key=0;
                twobit=0;
                j=i;
                while(j < end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOint[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j++;
                }
                if (insert_entry){
                    //Check for Ns
                    end = j;
                    for (j=end;j<end+4;j++)
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                    if (insert_entry){
                        currentPos[key]++;
                        total_num++;
                    }
                }
            }
            i++;
        }
    }

    //move file pointer to the beggining
    ifile.clear();
    ifile.seekg(0, ios::beg);

    //Reverse
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 22){
            if(chrom[i] == 'C' && chrom[i+1] == RevPAM){
                insert_entry = true;
                //Make the key
                end = i + 6;
                key=0;
                twobit=0;
                j=i+22;
                while(j > end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOintREV[chrom[j]-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j--;
                }
                //Check for Ns
                if (insert_entry){
                    end = j;
                    for (j=end;j>end-4;j--)
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                    if (insert_entry){
                        currentPos[key]++;
                        total_num++;
                    }
                }
            }
            i++;
        }
    }
    ifile.close();

    return total_num;
}

/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/******************************************   RESULT FUNCTIONS   *******************************************/
/***********************************************************************************************************/
/* Those functions populate the result table with the actual results. They get in by PAM and first the     */
/* forward strand and then the reverse strand.                                                             */
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/

/***********************************************************************************************************/
/* This function populates table B. Table A is already populated by the results of the count functions.    */
/* The input is the filename of the genome file, the PAM (NGG or NAG), a table with the positions that are */
/* currently occupied, table A, table B (to be filled), and a number (0-3).                                */
/* Essentially the function looks for hits exactly like the count functions above, but everytime one is    */
/* found instead of increasing a counter, the first 16 bases of the 20mer are looked up on table A and then*/
/* the result is inserted on the position indicated by table_A[16mer] + currentPos[16mer]. This temporary  */
/* table is identical with table A in structure but instead of the total number of hits it holds the number*/
/* of the hits of that specific 16mer that we already had inserted on table B. Hence, we know the exact    */
/* position to place the hit. The number is an ID that's not hardcoded to make the addition/removal of PAMS*/
/* simpler. It is added to each entry and used to identify the PAM when presenting the results to the user.*/
/***********************************************************************************************************/
void resultsNNG(const char *filename, char PAM, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID){
    unsigned int i,j,end,size,key,cnum;
    char RevPAM;
    string chrom(MAX_CHROM_LEN,'\0');
    string chromosome;
    MYTYPE twobit,value,rest4;
    bool insert_entry;

    //Depending on the input PAM find reverse for 2nd N
    if (PAM=='G')
        RevPAM = 'C';
    else if (PAM == 'A')
        RevPAM = 'T';
    else
        return;

    //Open the input file for reading
    ifstream ifile (filename,std::ifstream::in);

    //Add actual entries to the tables
    //Populating the hash table with forward entries first
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        //Read Header
        if (!isdigit(chrom[2]))
            chromosome.assign(&chrom[1],1);
        else
            chromosome.assign(&chrom[1],2);
        if(chromosome == "X") cnum = 23;
        else if (chromosome == "Y") cnum = 24;
        else if (chromosome == "M") cnum = 25;
        else cnum = atoi(chromosome.c_str());
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 22){
            if(chrom[i+21] == PAM && chrom[i+22] == 'G'){
                insert_entry = true;
                //Make the key
                end = i + 16;
                key=0;
                twobit=0;
                j=i;
                while(j < end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOint[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j++;
                }
                if (insert_entry){
                    //Make the value
                    twobit = 0;
                    value=0;
                    rest4 =0;
                    //Fix PAMID
                    value = PAMID;
                    //Fix strand (shift one - already 0)
                    value = value << STRANDBITS;
                    //Fix last nt before GG
                    value = ((value << LASTNTBITS)  | cTOint[(chrom[i+20])-ASCIIA]) & mask;
                    //Fix chromosome
                    value = ((value << CHROMBITS) | cnum) & mask;
                    //Fix position
                    value = ((value << POSBITS) | i) & mask;
                    //Fix last4 nt
                    end = j;
                    for (j=end;j<end+4;j++){
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                        twobit = cTOint[(chrom[j])-ASCIIA];
                        rest4 = ((rest4 << 2) | twobit) & mask;
                    }
                    if (insert_entry){
                        value = ((value << EXTENSIONBITS) | rest4) & mask;
                        //Add newly constructed pair to the structure
                        tableValues[hashTable[key] + currentPos[key]] = value;
                        currentPos[key]++;
                    }
                }
            }
            i++;
        }
    }

    //move file pointer to the beggining
    ifile.clear();
    ifile.seekg(0, ios::beg);

    //for (c=0;c<MAX_CHROM_NUM;c++){
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        if (!isdigit(chrom[2]))
            chromosome.assign(&chrom[1],1);
        else
            chromosome.assign(&chrom[1],2);
        if(chromosome == "X") cnum = 23;
        else if (chromosome == "Y") cnum = 24;
        else if (chromosome == "M") cnum = 25;
        else cnum = atoi(chromosome.c_str());
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 22){
            if(chrom[i] == 'C' && chrom[i+1] == RevPAM){
                insert_entry = true;
                //Make the key
                end = i + 6;
                key=0;
                twobit=0;
                j=i+22;
                while(j > end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOintREV[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j--;
                }
                if(insert_entry){
                    //Make the value
                    twobit = 0;
                    rest4 =0;
                    //Fix PAMID
                    value = PAMID;
                    //Fix strand
                    value = ((value << STRANDBITS) | 1) & mask;
                    //Fix last nt before GG
                    value = ((value << LASTNTBITS) | cTOintREV[(chrom[i+2])-ASCIIA]) & mask;
                    //Fix chromosome
                    value = ((value << CHROMBITS) | cnum) & mask;
                    //Fix position
                    value = ((value << POSBITS) | i) & mask;
                    //Fix last4 nt
                    end = j;
                    for (j=end;j>end-4;j--){
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                        twobit = cTOintREV[(chrom[j])-ASCIIA];
                        rest4 = ((rest4 << 2) | twobit) & mask;
                    }
                    if(insert_entry){
                        value = ((value << EXTENSIONBITS) | rest4) & mask;
                        //count newly constructed pair to the structure
                        tableValues[hashTable[key] + currentPos[key]] = value;
                        currentPos[key]++;
                    }
                }
            }
            i++;
        }
    }
    ifile.close();
}

/***********************************************************************************************************/
/* The input of this program must be a genome file which should have one line per chromosome that is       */
/* preceeded by a header line that contains only or starts with the chromosome name. The output is 2 files */
/* that are used by load memory to load all hits of all gRNAs of all PAMs into memory. The structure can be*/
/* queried efficiently by Results.o once loaded on memory.                                                 */
/***********************************************************************************************************/
int main (int argc, char *argv[]){
    char *usage = (char *) "\n\nUsage: %s -i <genome_filename> -o <output_full_path>\n\n";

    //Parse input
    int a;
    char *genome_file = NULL;
    char *out_file_path = NULL;
    while ((a = getopt (argc, argv, "i:o:")) != -1){
        if(a == 'i'){
            genome_file = optarg;
        }
        else if (a == 'o'){
            out_file_path = optarg;
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
            exit(1);
        }
    }

    if (genome_file == NULL || out_file_path==NULL) {
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        exit(1);
    }


    //Allocate space for the counters
    unsigned int *currentPos = (unsigned int*) calloc (KEY_NUM,sizeof(unsigned int));

    //Count all the results and the results per 16mer
    cout << "\nCounting NGG results. Please wait." << endl;
    unsigned int total_len=0;
    total_len = countNNG(genome_file,'G', currentPos);
    // cout << "Creating NGG index." << endl;

    //Allocate space equal to the total number of possible 16mers for the index table (table A)
    unsigned int *hashTable = (unsigned int*) calloc (KEY_NUM,sizeof(unsigned int));

    //Fill the index table A with the counters per 16mer
    hashTable[0]=0;
    unsigned long v;
    for(v=1;v<KEY_NUM;v++){
        hashTable[v] = hashTable[v-1] + currentPos[v-1];
    }

    //Allocate space for the result table
    MYTYPE *tableValues = new MYTYPE [total_len];

    //re allocating the counters (so that they'll be 0's). This will keep the
    //number of results put in table B at every step of the populating functions
    free(currentPos);
    currentPos = (unsigned int*) calloc (KEY_NUM,sizeof(unsigned int));

    //Fill table B with the results
    cout << "\nCalculating NGG results. Please wait." << endl;
    resultsNNG(genome_file,'G',currentPos,tableValues,hashTable,0);
    // cout << "Calculating NGG count." << endl;

    //Write index and result tables on binary files
    cout << "\nWriting tables on index.bin and data.bin" << endl;
    FILE *f1,*f2;
    char index_file[strlen(out_file_path) + strlen("/index.bin") + 1], data_file[strlen(out_file_path) + strlen("/data.bin") + 1];
    strcpy(index_file,out_file_path);
    strcpy(data_file,out_file_path);
    if (out_file_path[strlen(out_file_path) - 1] == '/'){
        f1 = fopen ( strcat(index_file,"index.bin"), "wb" );
        f2 = fopen ( strcat(data_file,"data.bin"), "wb" );
    }
    else{
        f1 = fopen ( strcat(index_file,"/index.bin"), "wb" );
        f2 = fopen ( strcat(data_file,"/data.bin"), "wb" );
    }
    fwrite(hashTable,sizeof(unsigned int),KEY_NUM,f1);
    fwrite(tableValues,sizeof(MYTYPE),total_len,f2);
    fclose(f1);
    fclose(f2);

    cout << "\nTables created.\n" << endl;

    return 1;
}

