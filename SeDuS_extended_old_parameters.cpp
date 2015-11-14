/****************************************************************************
 **                                                                        **
 **  SEDUS, Segmental Duplication Simulator                                **
 **  Copyright (C) 2015 Diego A. Hartasánchez, Marina Brasó-Vives,         **
 **  Oriol Vallès-Codina, Juanma Fuentes-Díaz and Arcadi Navarro,          **
 **  Institut de Biologia Evolutiva UPF-CSIC                               **
 **                                                                        **
 **  This file is part of SEDUS.                                           **
 **                                                                        **
 **  SEDUS is free software: you can redistribute it and/or modify         **
 **  it under the terms of the GNU General Public License as published by  **
 **  the Free Software Foundation, either version 3 of the License, or     **
 **  (at your option) any later version.                                   **
 **                                                                        **
 **  SEDUS is distributed in the hope that it will be useful,              **
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of        **
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
 **  GNU General Public License for more details.                          **
 **                                                                        **
 **  You should have received a copy of the GNU General Public License     **
 **  along with this program.  If not, see http://www.gnu.org/licenses/.   **
 **                                                                        **
 ****************************************************************************
 **          Authors: Diego A. Hartasánchez, Marina Brasó-Vives,           **
 **                   Juanma Fuentes-Díaz, Oriol Vallès-Codina,            **
 **                   and Arcadi Navarro                                   **
 **  Website/Contact: http://www.biologiaevolutiva.org/sedus/              **
 **             Date: Jul 28 2015                                          **
 **          Version: 1.10                                                 **
 ****************************************************************************/

/* This is the command-line version of SeDuS */

/*ARGUMENTS: SimulationID */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <ctime>
#include <algorithm>

#include <cstddef>
#include <array>
#include <vector>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

using namespace std;

////////////////////////////////////////
//// SIMULATION PRINCIPAL VARIABLES ////
////////////////////////////////////////

int how_many_new_mutations_in_dup = 0; // WHC: counting how many new mutations happened in block_0, 2 or 4


int N = 100; // Population size

// WHC: PROMETHEUS should be an even number
int PROMETHEUS = 100; // Number of generations for each genealogy

int SUPERTIME = 2; // Number of simulations per execution
int BLOCKLENGTH = 10000; // Block length
int SAMPLE = 50; // Sample size
#define MUTTABLESIZE 1000000 // Maximum number of mutations (size of muttable)

// #define B 3 // Maximum number of blocks per chromosome
// #define B 5 // WHC: Maximum number of blocks per chromosome
#define B 6 // WHC: Maximum number of blocks per chromosome

// WHC: define index for each block, instead of using 0, 1, 2 for original, single-copy and duplicated
// int ori_index = 0, single_1 = 1, dup_1 = 2, single_copy_2 = 3, dup_2 = 4;
int ori_index = 0, single_1 = 1, spacer = 2, dup_1 = 3, single_copy_2 = 4, dup_2 = 5;

#define numOfBins 5 // Number of segments in which we divide each block (exclusively used for analysis)
// #define maxNumOfHS 10  // Maximum number of crossover hotspots
// #define maxNumOfHS 3  // Maximum number of crossover hotspots
#define maxNumOfHS 5  // Maximum number of crossover hotspots
#define numofsamples 1 // Number of different samples taken per run

// int BIGTIME = 70;  // BIGTIME * N is the total number of generations simulated per run
int BIGTIME = 70;  // WHC: BIGTIME * N is the total number of generations simulated per run, UNTIL phaseIV()
int BURNINTIME = 20; // BURNINTIME * N is the number of generations in Phase I
// int STRUCTUREDTIME = 20; // STRUCTUREDTIME * N is the number of generations in Phase II
int STRUCTUREDTIME = 30; // WHC: STRUCTUREDTIME * N is the number of generations in Phase II

// WHC: STRUCTURED_2_TIME * N is the number of generations in Phase IV
int STRUCTURED_2_TIME = 30;


// WHC: PHASE_V_TIME * N is the number of generations of Phase V
int PHASE_V_TIME = 30;
int PHASE_V_LENGTH = (int) PHASE_V_TIME * N;

// WHC: phaseVI time
int PHASE_VI_TIME = 60;
int PHASE_VI_LENGTH = (int) PHASE_VI_TIME * N;

int TIMELENGTH = (int) BIGTIME * N; // Total number of generations (including all phases) (not including phaseIV(), phaseV(), phaseVI())
int BURNIN = (int) BURNINTIME * N; // Number of generations in phase I
int STRUCTURED = (int) STRUCTUREDTIME * N; // Number of generations in phase II

// WHC: Number of generations in phase IV
int STRUCTURED_2 = (int) STRUCTURED_2_TIME * N;

float THETA = 0.001; // Population scaled mutation rate: THETA = 4 N mu
//float R = 10; // Population scaled crossover rate: R = 4 N rho*BLOCKLENGTH
float R = 100;
float C = 0.5; // Population scaled gene conversion rate: C = 4N*kappa*lambda

// int numHS = 1; // Number of hotspots
// int numHS = 3; // WHC: Number of hotspots, single_2 is a scaled hopspot representing a 2MB region
int numHS = 5; // WHC: added spacer_1, spacer_2; spacer_1 = block 2, spacer_2 = block 4; to eliminate the effects of balancing selection on single_copy

// WHC: hard cord to initiate these values
//int crossoverBegin[maxNumOfHS] = {0, 3 * BLOCKLENGTH, 4 * BLOCKLENGTH}; // Start point of crossover regions
//int crossoverEnd[maxNumOfHS] = {3 * BLOCKLENGTH, 4 * BLOCKLENGTH, 5 * BLOCKLENGTH}; // End point of crossover regions

int crossoverBegin[maxNumOfHS] = {0, 2 * BLOCKLENGTH, 3 * BLOCKLENGTH, 4 * BLOCKLENGTH, 5 * BLOCKLENGTH}; // Start point of crossover regions
int crossoverEnd[maxNumOfHS] = {2 * BLOCKLENGTH, 3 * BLOCKLENGTH, 4 * BLOCKLENGTH, 5 * BLOCKLENGTH, 6 * BLOCKLENGTH}; // End point of crossover regions

// WHC: artificially selected values, to represent 2MB region as crossover hotspot

// WHC: REMEMBER!!! need to change crossoverRatio[] in int main(), before phaseIV()
// WHC: not functional any more
float crossoverRatio[maxNumOfHS] = {0.08, 0.04, 0.04, 0.80, 0.04}; // Relative weights of crossover regions

string letter = ""; // Simulation ID
string str;

////////////////////////////////
//// DECLARATION OF CLASSES ////
////////////////////////////////

struct chrom { // Chromosomes
  int b; // Number of blocks
  int mpb[B]; // Number of mutations per block
  std::vector<std::vector<int>> mutation;
  chrom(){
    mutation.resize(B);
    for(unsigned int i=0;i<mutation.size();i++){mutation[i].resize(BLOCKLENGTH/2);} // WHC: is the size just for saving memory? YES
    //    for(unsigned int i=0;i<mutation.size();i++){mutation[i].resize(BLOCKLENGTH);} // WHC: is the size just for saving memory
  }
};

struct mutation { // Mutations
  int position;
  int block;
  float frequency;
};

struct prev_pres {
  int prev;
  int pres;
};

struct fertility_info{
  int x;
  int y;
  bool recombinatrix;		// WHC: useless?
};


///////////////////////////////
//// VARIABLES DECLARATION ////
///////////////////////////////

//// PARAMETERS OF EVENT ////
float mu = THETA / (4 * N); // Mutation rate per nucleotide and generation
float rho = R / (4 * N * BLOCKLENGTH); // Crossover rate per nucleotide and generation

// WHC: may consider changing
float meanTractLength = 100; // Mean Gene Conversion Tract Length (lambda)

float kappa = C / (4 * N * meanTractLength); // Gene Conversion Initiation Rate

// WHC: may consider changing
float meps = 0; // Length of the 100% identity tract
float similarityInConvTract = 0; // Percent of similarity required for conversion in all conversion tract

// float donorRatio = 0.5; // Percentage of gene conversion events that occur from the original to the duplicated block

/* WHC: donorRatio should be a 2-D array, donorRatio[ori][dup_1] = the percent of donor is ori compared with is dup_1, i.e.
 * (dis-enable command-line input of donorRatio for now
 * (should actually be 5x5 array, as ori_index = 0, dup_1 = 2, dup_2 = 4

 *                ori              single_1     dup_1         single_2         dup_2
 *     ori         -1                 -1      x(ori->dup_1)      -1         x(ori->dup_2)
 *     single_1    -1                 -1         -1              -1              -1
 *     dup_1     1 - x(ori->dup_1)    -1         -1              -1         x(dup_1->dup_2)
 *     single_2    -1                 -1         -1              -1              -1
 *     dup_2     1 - x(ori->dup_2)    -1    1 - x(dup_1->dup2)   -1              -1
 *
 *     correspondingly, there should be a step in void conversion(), to pick 2 duplications first
 *     i.e. if (p <= 1/3) choose ori & dup_1; else if (1/3 < p <= 2/3) choose ori & dup_2; else if (p > 2/3) choose dup1&2;
 */

const float ori_to_dup_1 = 0.5, ori_to_dup_2 = 0.5, dup_1_to_dup_2 = 0.5;
// Percentage of gene conversion events that occur from the original to the duplicated block
// WHC: donorRatio initialization
float donorRatio[5][5] = {
  {-1, -1, ori_to_dup_1, -1, ori_to_dup_2},
  {-1, -1, -1, -1, -1},
  {1 - ori_to_dup_1, -1, -1, -1, dup_1_to_dup_2},
  {-1, -1, -1, -1, -1},
  {1 - ori_to_dup_2, -1, 1 - dup_1_to_dup_2, -1, -1}
};




// WHC: this process of choose same/different chroms could remain the same
// WHC: sameDifIGC, when == 1, mean all from same chrom; if == 0, means all from different chrom
float sameDifIGC = 1; // Proportion of IGC events that occur between copies in the same chromosome (1- between copies in homologous chromosomes)

bool dupType = 1; // Duplication mechanism: 0 = to the same chromosome / 1 = to the partner chromosome

//// GENEALOGICAL MATRICES ////
std::vector<std::vector<int>> ancestry; // Matrix codifying the ancestor of each chromosome in each generation
std::vector<fertility_info> fertility_list;
std::vector<fertility_info> fertility_list_ini;
std::vector<std::vector<bool>> fertility;// Boolean Matrix codifying if a chromosome has descendants in the final generation of each era (each era is equivalent to PROMETHEUS generations) or not
std::vector<std::vector<bool>> fertility_ini;//Matrix fertility initialization
std::vector<std::vector<bool>> recombimatrix;// Boolean Matrix codifying if a chromosome comes from a recent recombination or not (has one or two parents)
std::vector<std::vector<int>> duplicontent;// Array indicating which chr carry the duplication (with present and previous lines)

// WHC: duplicontent for phaseIV()
vector<vector<int>> duplicontent_2;

// WHC: lose_of_duplicontent[] for phaseVI()
vector<vector<int>> lose_of_duplicontent;

std::vector<std::vector<bool>> IGCmatrix;// Boolean Matrix codifying if a chromosome undergoes IGC or not 

//// GLOBAL VARIABLES AND STATISTICAL QUANTITIES ////
std::vector<chrom> table;
std::vector<std::vector<chrom*>> pointer; // GENOMIC INFORMATION
int era, run; // t = generation inside each era, era = PROMETHEUS generations inside each phase (from 0 to BIGTIME-1), run = number of simulation runs (from 0 to SUPERTIME-1)

// WHC: maybe too short for ours
// int fixationTrajectory[20000 + 1]; // Absolute frequency of the duplication in each generation of fixation process
// WHC: only allow 30 * 1000?? 30->structured, 1000->#individuals
int fixationTrajectory[30000 + 1]; // Absolute frequency of the duplication in each generation of fixation process



// WHC: remember to change before use
int phaseVI_trajectory[30000 + 1];

std::vector<bool> multihit; // Record the positions in which a mutation has occurred

// WHC: multihit for single-copy block
vector<bool> multihit_single;

int duplicationFreq; // Absolute frequency of the duplication in the present generation
bool duFreq; // Duplication has occurred or not

// WHC: Duplication_2 has occurred or not
bool duFreq_2;

// WHC: whether or not loosing dup_1 has occured
bool loseFreq;

//// GLOBAL INTERNAL VARIABLES FOR FSL ////
int MutCount; // Total number of mutations segregating (the fixed ones are erased) in each moment
struct mutation muttable[MUTTABLESIZE], temporalmuttable[MUTTABLESIZE]; // Register of all the mutations

//// SAMPLING ////
std::vector<int> sample; // Randomly sampled individuals
int sampleN[] = {SAMPLE}; // Size of the samples (only one sample in this case)

//// FILES ////
ofstream profile;// N PROMETHEUS SUPERTIME BLOCKLENGTH SAMPLE MUTTABLESIZE B BIGTIME BURNIN/N STRUCTURED/N mu rho kappa meanTractLength
//ofstream auxx; // Auxiliary print file that can be used for particular needs (open and close lines must be uncommented)
ofstream samplefile[B + 2][numofsamples][2]; // For each block and the collapsed. 0 = pi; 1 = S. For each sample
ofstream mutationsNewFile[B+1]; // new mutation file, with ms-like format
ofstream SFS[B+2]; // site frequency spectra for each block + collapsed only for last era

// WHC: ofstream for inter-block pairwise divergence
ofstream interblock_divergence[3]; // 0 -> between original and dup_1, 1 -> between original and dup_2, 2 -> between dup_1 and dup_2

///////////////////
//// FUNCTIONS ////
///////////////////

struct prev_pres phaseI();
void phaseII(int,int,int, float);
void phaseIII(float);

// WHC: for phaseIV(); it is a similar function to phaseII()
void phaseIV(int, int, int, float);

// WHC: for phaseV(); it is a similar function to phaseIII()
void phaseV(float);

// WHC: for phaseVI(); similar to phaseIV(), but is similating losing dup_1
void phaseVI(int, int, float);

void open_files(); // Opening files
void close_files(); // Closing files

void genealogy(float, int, float); // rho, 0/1(non structured or structured) //// Filling genealogy matrices in each PROMETHEUS
void genealogy_2(float, int, float);
void parentpicking(int[maxNumOfHS], int[maxNumOfHS], float[maxNumOfHS], int, int, int,int,int); // crossoverBegin, crossoverEnd //// Create new generation from previous one (with recombination)

// WHC: parentpicking for phaseVI
void parentpicking_for_phaseVI(int[maxNumOfHS], int[maxNumOfHS], float[maxNumOfHS], int, int, int,int,int);

void duplication(int,int,bool); // Create Duplication for eva (first duplicated chromosome)

// WHC: copy of duplication() for phaseIV()
void duplication_2(int,int,bool); // Create Duplication for eva (first duplicated chromosome)

// WHC: losing dup_1 in phaseVI
void lose_of_duplication(int, int);

void mutation(float, int, int); // For each fertile chromosome decide if a mutation happens and execute it if necessary

// WHC: mutation() for phaseVI
void mutation_for_phaseVI(float, int, int);

// WHC: will use pointer donorRatio instead of float
void conversion(float, int, int, int, float (*donorRatio)[5], float); // Only when duFreq is true. For each fertile chromosome decide if conversion happens and execute it if necessary
void conversion_for_phaseVI(float, int, int, int, float (*donorRatio)[5], float);

// WHC: conversion() for phaseVI()

void statistics(int, bool); // Execution of all the statistic calculations (for each Era)

// WHC: statistic() for phaseVI
void statistics_for_phaseVI(int, bool);

void FSL(int); // Count of All the Segregating, Fixed, Lost and Shared sites

// WHC: a test FSL(), to solve the problem that multihit[] == true just grows up since phaseIV()
void FSL_since_phaseIV(int);

void copychr(int, int, int, int); // Copy a chromosome (when there is no recombination)

// WHC: new copychr() for phaseVI
void copychr_for_phaseVI(int, int, int, int);

int location(int, int, int, int); // Returns the location of a mutation or point in the mutation vector of a chromosome
void EraseFixedMutations(int, int, int); // Erase fixed mutations from chromosomes that have them (not necessary to consider them any more)
int SearchMutation(int, int, int); // Search a given mutation inside muttable (returns its position)

float * SiteFrequencySpectrumPrint(int, int, int, bool); // Calculate pi and S values for a block

// WHC: SiteFrequencySpectrumPrint() for phaseVI
float * SiteFrequencySpectrumPrint_for_phaseVI(int, int, int, bool);


void DivergenceForAll(int, int, int); // Average divergent positions between two blocks in the same chromosome

int DupliFreq(int, int, int); // Calculate Duplication frequency in the population

// WHC: DupliFreq() for phaseVI
int DupliFreq_for_phaseVI(int, int, int);

int muFrequencyIntWholePopAndSample(int, int, int, int); // Calculate absolute frequency of a given mutation
int muFrequencyCollapsedCallingFromSample(int, int, int); // Calculates absolute frequency of a given mutation collapsing both blocks from each individual

int tractpql(float); // Returns a given random tract length from the mean tract length
void SamplingIndividuals(int); // Sample corresponding number of individuals (register their number in sample[])
int GenerateFixationTrajectory(int, int); // Generates fixation trajectory of the duplication (fixed or not)

// WHC: for generating trajectory for phase VI
void Generate_phaseVI_trajectory(int);

// WHC: randomly pick 0, 1, 2 with desired propotions
int pick_a_pair(int);


void print_fertility();
float round(float, int);
int minim (int, int);
bool sortx (fertility_info,fertility_info);
bool sorty (fertility_info,fertility_info);

float number_of_nucleotide_diff(int, int, int, int, int);
float pairwise_divergence_between(int, int, int, int);

//////////////
//// MAIN ////
//////////////

int main ( int argc, char* argv[] ) { // WHC: argc is the # of arguments passed from command line

  //  cout << sizeof(phaseVI_trajectory) / sizeof(phaseVI_trajectory[0]) << endl;
  //  cout << sizeof(fixationTrajectory) / sizeof(fixationTrajectory[0]) << endl;
  if (sizeof(phaseVI_trajectory) / sizeof(phaseVI_trajectory[0]) < PHASE_VI_LENGTH) {
    cout << "phaseVI_trajectory[] length too small, need to be longer!\n";
    exit(0);
  }

  if (sizeof(fixationTrajectory) / sizeof(fixationTrajectory[0]) < STRUCTUREDTIME || sizeof(fixationTrajectory) / sizeof(fixationTrajectory[0]) < STRUCTURED_2_TIME) {
    cout << "fixationTrajectory[] length too small.\n";
    exit(0);
  }
  using namespace std;
  int correctArguments = 0;
   
  int argN, argk, argBlockLength, argSample, argBurnIn, argBigTime;
  int sup, timeToFixation=0, num;  
  float argTheta, argR, argC, argMTL, argMEPS, argDonorRatio, argSameDifIGC, Beg, End, HSRatio;
  //  bool notSC = 0;
  bool notSC = 1;
		
  if (argc < 2) { // Check the value of argc. If not enough parameters have been passed, inform user and exit.
    cout << "You need to provide a Simulation ID\n";
    exit(0);
  } else { // IF AT LEAST SimulationID IS PROVIDED
    letter = argv[1]; // WHC: letter is SimuID
    for (int i = 2; i < argc; i++) { // ITERATE OVER THE REST OF ARGUMENTS
			
      if (i + 1 != argc){ // Check that it hasn't finished parsing
	std::string str = argv[i];
	const char *cstr = str.c_str(); // WHC: return a C-type string pointer
	if (strcmp(cstr,"-s") == 0) { //SUPERTIME
	  sscanf (argv[i+1],"%d",&sup); SUPERTIME = sup; // WHC: transform argv[i+1] to form %d and save into int sup
	} else if (strcmp(cstr,"-n") == 0) {// POPULATION SIZE (N)
	  sscanf (argv[i+1],"%d",&argN); N = argN;
	} else if (strcmp(cstr,"-k") == 0) {// PROMETHEUS (GENERATIONS BETWEEN SNAPSHOTS K)
	  sscanf (argv[i+1],"%d",&argk); PROMETHEUS = argk;
	} else if (strcmp(cstr,"-b") == 0) {// BLOCKLENGTH
	  sscanf (argv[i+1],"%d",&argBlockLength); BLOCKLENGTH = argBlockLength;
	} else if (strcmp(cstr,"-z") == 0) {// SAMPLE SIZE
	  sscanf (argv[i+1],"%d",&argSample); SAMPLE = argSample;
	} else if (strcmp(cstr,"-i") == 0) {// BURNIN
	  sscanf (argv[i+1],"%d",&argBurnIn); BURNINTIME = argBurnIn;
	} else if (strcmp(cstr,"-t") == 0) {// BIGTIME
	  sscanf (argv[i+1],"%d",&argBigTime); BIGTIME = argBigTime;
	} else if (strcmp(cstr,"-p") == 0) {// 	DUPLICATION TYPE (SAME / HOMOLOGOUS CHROMOSOME)
	  dupType = 0;
	} else if (strcmp(cstr,"-u") == 0) {// MUTATION RATE THETA
	  sscanf (argv[i+1],"%f",&argTheta); THETA = argTheta;
	} else if (strcmp(cstr,"-r") == 0) {//CROSSOVER RATE R
	  sscanf (argv[i+1],"%f",&argR); R = argR;
	} else if (strcmp(cstr,"-c") == 0) { //IGC RATE C
	  sscanf (argv[i+1],"%f",&argC); C = argC;
	} else if (strcmp(cstr,"-l") == 0) { // MEAN TRACT LENGTH
	  sscanf(argv[i+1], "%f", &argMTL); meanTractLength = argMTL;
	} else if (strcmp(cstr,"-m") == 0) { // MINUMUM EFFICIENT PROCESSING SEGMENT
	  sscanf(argv[i+1], "%f", &argMEPS);  meps = argMEPS;
	} else if (strcmp(cstr,"-y") == 0) { // RATIO IGC SAME / HOMOLOGOUS CHROMOSOMES
	  sscanf(argv[i+1], "%f", &argSameDifIGC);  sameDifIGC = argSameDifIGC;
	} else if (strcmp(cstr,"-d") == 0) { // DONOR RATIO
	  // WHC: dis-enable for now
	  cout << "You must hard-code donor ratio matrix now!\n";
	  exit(0);
	  //	  sscanf (argv[i+1], "%f", &argDonorRatio); donorRatio = argDonorRatio;
	  //	  if(donorRatio < 0 || donorRatio > 1 ) { cout << "Donor ratio must be between 0 and 1.\n"; exit(0);}
	} else if (strcmp(cstr,"-f") == 0) { // TIME TO FIXATION
	  sscanf (argv[i+1],"%d",&timeToFixation);
	} else if (strcmp(cstr,"-w") == 0) { // WHOLE-REGION CROSSOVER
	  crossoverBegin[0] = 1; crossoverEnd[0] = B*BLOCKLENGTH; crossoverRatio[0] = 1; notSC = 1;
	} else if (strcmp(cstr,"-h") == 0) { // HOTSPOT CROSSOVER
	  // WHC: disable this function for now
	  cout << "Sorry but this function is disabled now.\n";
	  exit(0);
	  sscanf (argv[i+1],"%d",&num); numHS = num; notSC = 1;
	  if (argc < (i+4+(num-1)*3)) { cout << "You need to provide numHS and HS_start HS_end HS_ratio for each hotspot.\n"; exit(0);}
	  float sumRatio = 0.000;
	  for (int HS = 0; HS < numHS; HS++) {
	    sscanf (argv[i+2+HS*3],"%f",&Beg); crossoverBegin[HS] = (int) (Beg * BLOCKLENGTH); // WHC: read crossover hopspots Begin & End positions into 2 array
	    sscanf (argv[i+3+HS*3],"%f",&End); crossoverEnd[HS] = (int) (End * BLOCKLENGTH);
	    if(Beg < 0 || End > 3 || Beg > End) {
	      cout << "Crossover regions should be within the range [0," << B << "] and HS_start should be smaller than HS_end.\n"; exit(0);
	    }
	    sscanf (argv[i+4+HS*3],"%f",&HSRatio); crossoverRatio[HS] = HSRatio;
	    sumRatio = sumRatio + HSRatio;
	  }
	  if(sumRatio != 1.000) { cout << "Sum of hotspot ratios must be 1.\n"; exit(0);}
				 
	}
      }
    }
  }

  
  multihit.resize(BLOCKLENGTH); // WHC: allocate BLOCKLENGTH-long size for vector multihit
  multihit_single.resize(BLOCKLENGTH);
  
  table.resize(4*N);	      // WHC: a table of struct chroms
  pointer.resize(2);for(unsigned int i=0;i<pointer.size();i++){pointer[i].resize(2*N);} // WHC: pointers to chroms
  sampleN[0] = {SAMPLE};
  sample.resize(2*N);	// WHC: individuals
  ancestry.resize(2 * N);for(unsigned int i=0;i<ancestry.size();i++){ancestry[i].resize(PROMETHEUS);}
  fertility.resize(2 * N);for(unsigned int i=0;i<fertility.size();i++){fertility[i].resize(PROMETHEUS);}
  fertility_ini.resize(2 * N);for(unsigned int i=0;i<fertility_ini.size();i++){fertility_ini[i].resize(PROMETHEUS);}
  recombimatrix.resize(2 * N);for(unsigned int i=0;i<recombimatrix.size();i++){recombimatrix[i].resize(PROMETHEUS);}
  duplicontent.resize(2);for(unsigned int i=0;i<duplicontent.size();i++){duplicontent[i].resize(2*N);}

  duplicontent_2.resize(2);for(unsigned int i=0;i<duplicontent_2.size();i++){duplicontent_2[i].resize(2*N);}

  // WHC: lose_of_duplicontent[] records which chrom only has 4 blocks instead of 5
  // (chrom at the beginning carrying 4 blocks, similar to duplicontent[]
  lose_of_duplicontent.resize(2); for (unsigned int i = 0; i < lose_of_duplicontent.size(); ++i) {lose_of_duplicontent[i].resize(2 * N);}
  
  IGCmatrix.resize(2 * N);for(unsigned int i=0;i<IGCmatrix.size();i++){IGCmatrix[i].resize(PROMETHEUS);}

  correctArguments = 1;
  //  if(timeToFixation > 20*N) { cout << "Time to fixation must be smaller than 20N.\n"; exit(0);}
  if(timeToFixation > 50*N) { cout << "Time to fixation must be smaller than 50N.\n"; exit(0);}
  mu = THETA / (4 * N);
  rho = R / (4 * N * BLOCKLENGTH);
  kappa = C / (4 * N * meanTractLength);

  //  kappa = 0;
  
  TIMELENGTH = (int) BIGTIME * N; 
  BURNIN = (int) BURNINTIME * N; 


  //// SET DEFAULT PARAMETERS
  // WHC: not default notSC == 1
  if(notSC == 0){		// WHC: mode is SCC, crossover happens only in single-copy region
    cout << "notSC should be 1\n";
    exit(0);
    crossoverBegin[0] = BLOCKLENGTH;
    crossoverEnd[0] = 2*BLOCKLENGTH;
    crossoverRatio[0] = 1;
  }
	
  if (correctArguments == 1){ // WHC: to set random number generator from current TIME
    int i, j, h, o; //endTime;
    time_t seconds;	// WHC: seconds between now and 1970/1/1
    time(&seconds);	// WHC: get the abovementioned seconds
    srand((unsigned int) seconds);

    open_files();	// WHC: set up output file streams

    //PROFILE IN TAB FORMAT
    profile << "Runs\tSampleSize\tGenerationsBetweenSnapshots(k)\tPopulationSize(N)\tTheta\tBlockLength\tDuplicationOrigin\tBurnIn\tTimeToFixation\t";
    profile  << "TimeLength\tC\tMeanTractLength\tDonorRatio\tMEPS\tProportionSameDifIGC(w)\tR\tNumOfHS";
    for (j = 0; j < numHS; j++) {
      profile  << "\tHS" << j << "_st\tHS" << j << "_end\tHS" << j <<"_ratio";
    } 
    profile << "\n";
		
    profile << SUPERTIME << "\t" << SAMPLE << "\t" << PROMETHEUS << "\t" << N << "\t" << THETA << "\t" << BLOCKLENGTH << "\t" << dupType << "\t" ;
    profile << BURNIN << "\t" << timeToFixation << "\t" << TIMELENGTH << "\t" << C << "\t"<< meanTractLength << "\t";
    profile << donorRatio << "\t" << meps  << "\t" << sameDifIGC << "\t" << R << "\t" << numHS << "";
    for (j = 0; j < numHS; j++) {
      profile << "\t" << crossoverBegin[j] << "\t" << crossoverEnd[j] << "\t" << crossoverRatio[j] << "";
    }
    profile << "\n";

    //  for (j = 0; j < B; j++) {
    //  	mutationsNewFile[j] << "ms " << SAMPLE << " " << SUPERTIME << " -s 5\n" << seconds << "\n";
    //  }

    //fertility initialization
    for (int i=0 ; i < 2*N ; i++){ // WHC: Except for the last generation, all fertility = false
      fertility_ini[i][PROMETHEUS-1] = true;
      for (int tt=0 ; tt < PROMETHEUS-1 ; tt++) {fertility_ini[i][tt] = false;}
    }

    
    // INITIALIZATION SUPERTIME
    for (run = 0; run < SUPERTIME; run++) {
      cout << "SUPERTIME = " << run << "\n";
      time(&seconds); // WHC: second time, possibly redundency
      //srand((unsigned int) seconds);
      srand((unsigned int) seconds);
			

      ////////////////////////////////////
      //  BUILD THE INITIAL POPULATION  //
      ////////////////////////////////////
      for (i = 0; i < 4 * N; i++) {
	table[i].b = 2; // Each chromosome begins with 2 blocks (b)
	for (j = 0; j < B; j++) {
	  table[i].mpb[j] = 0; // For each block Mutations per block = 0
	}
      }
      for (h = 0; h < 2; h++) {
	for (i = 0; i < 2 * N; i++) {
	  pointer[h][i] = &table[i + 2 * N * h]; // Inside pointer <- table values
	}
      }

      MutCount = 0;
      duFreq = false;

      duFreq_2 = false;

      loseFreq = false;
      
      for (j = 0; j < BLOCKLENGTH; j++) { multihit[j] = false; multihit_single[j] = false; }

      //////////////////////
      //////// RUN /////////
      //////////////////////

// WHC: hard cord to initiate these values
//      crossoverBegin[maxNumOfHS] = {0, 3 * BLOCKLENGTH, 4 * BLOCKLENGTH}; // Start point of crossover regions
//      crossoverEnd[maxNumOfHS] = {3 * BLOCKLENGTH, 4 * BLOCKLENGTH, 5 * BLOCKLENGTH}; // End point of crossover regions
      // WHC: before phaseIV(), crossover could only happen in block 0, 1, 2
      //      crossoverRatio[maxNumOfHS] = {1, 0, 0}; // Relative weights of crossover regions
      crossoverRatio[0] = 1;
      crossoverRatio[1] = 0;
      crossoverRatio[2] = 0;
      crossoverRatio[3] = 0;
      crossoverRatio[4] = 0;
      
      /*  PHASE I: BURN-IN  */
      cout << "PHASE I" << "\n";
      prev_pres ret = phaseI();
      /* END PHASE I */


      crossoverRatio[0] = 0.50;
      crossoverRatio[1] = 0.25;
      crossoverRatio[2] = 0.25;
      crossoverRatio[3] = 0;
      crossoverRatio[4] = 0;

      
      /*  PHASE II: STRUCTURED TRAJECTORY  */
      cout << "PHASE II" << "\n";
      //endTime=
      phaseII(timeToFixation,ret.prev,ret.pres, kappa);
      /* END PHASE II */

      /*  PHASE III: RESOLUTION  */
      cout << "PHASE III" << "\n";
      phaseIII(kappa);
      /* END PHASE III */


      //      crossoverBegin[maxNumOfHS] = {0, 3 * BLOCKLENGTH, 4 * BLOCKLENGTH}; // Start point of crossover regions
      //      crossoverEnd[maxNumOfHS] = {3 * BLOCKLENGTH, 4 * BLOCKLENGTH, 5 * BLOCKLENGTH}; // End point of crossover regions

      //      crossoverRatio[maxNumOfHS] = {0.15, 0.8, 0.05}; // Relative weights of crossover regions

      // need to change this here
      
      crossoverRatio[0] = 0.08;
      crossoverRatio[1] = 0.04;
      crossoverRatio[2] = 0.04;
      crossoverRatio[3] = 0.80;
      crossoverRatio[4] = 0.04;
      
      /* PHASE IV: STRUCTURED_2 TRAJECTORY */
      // WHC: for generating dup_2
      cout << "PHASE IV" << '\n';
      // WHC: I do not think ret.prev, ret.pres is necessary...
      phaseIV(timeToFixation,0, 1, kappa);
      /* END PHASE IV */

      /* PHASE V: RESOLUTION after 2nd duplication */
      cout << "PHASE V" << '\n';
      phaseV(kappa);

      /* PHASE VI: LOSING dup_1 */
      cout << "PHASE VI" << '\n';
      phaseVI(0, 1, kappa);
      
      for (j = 0; j < B+2; j++) {
      	for (o = 0; o < numofsamples; o++) {samplefile[j][o][0] << "\n";samplefile[j][o][1] << "\n";}
      }
      for (int k = 0; k < 3; ++k) {
	interblock_divergence[k] << '\n';
      }
      cout << "Finished Run" << "\n";
      // auxx << endTime <<"\n";

    } // END SUPERTIME
    close_files();
  }
  else { // INCORRECT ARGUMENTS
    cout << "INCORRECT ARGUMENTS\nUsage: SimulationID R C F/NF SC/WR/HS\nF/NF: Fixed or not fixed duplication fixation trajectory\nIn case of HS you should include the number of hostspots and each crossover begin point and crossover end point as arguments\n";
  }
  // clock_t end = clock();
  //  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  // cout << "time: " << elapsed_secs;

  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////  FUNCTIONS   /////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

struct prev_pres phaseI(){
  // 0 -> BURNIN/PROMETHEUS

  // WHC: I think prev and pres is for switching between present generation and previous generations within pointer[][]
  // WHC: so pointer will always maintain 2 generations' information at a time, in parallel, and switch pres & prev then
  
  int prev = 0;
  int pres = 1;
  bool does_print = false;
  prev_pres ret;
  for (era = 0; era < (int) BURNIN / PROMETHEUS; era++) {
    genealogy(rho * BLOCKLENGTH, 0, 0);// GENEALOGY (filling recombimatrix, ancestry and fertility matrices, determines all population genealogy)
    int prom=-1;
    prev = 0;
    pres = 1;
    for (std::vector<fertility_info>::iterator it=fertility_list.begin(); it!=fertility_list.end(); ++it){
      int i = (*it).x;
      int t = (*it).y;
      if(prom!=t){
	if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
      }
      prom = t;
      // WHC: parentpicking will finish recombinations between father & partner, and save results to pointer[pres][i]
      parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS, prev,pres,i,t);
      mutation(mu, i, pres);
    }
    // CALCULATE THE STATISTICS
    statistics(pres, does_print);
  }
  ret.prev = prev;
  ret.pres = pres;
  return ret;
}

void phaseII(int timeToFixation,int prev, int pres, float k){
  bool does_print = false;
  // Generating Fixation Trajectory in a diploid population
  int endTime = GenerateFixationTrajectory(STRUCTURED + 1, timeToFixation);
  // Picking the first chromosome with the duplication (eva)
  int eva = (int) (rand() % (2 * N));
  for (int i = 0; i < 2 * N; i++) {
    duplicontent[0][i] = i;
    duplicontent[1][i] = i;
  }
  duplicontent[1][0] = eva;
  duplicontent[1][eva] = 0;
  // Duplicate eva (create a new block with old mutations...)
  duplication(eva,pres,dupType); // WHC: the pres (present) here, will be used as prev (previous) in duplication(); because this is the start of PhaseII
  duFreq = true;

  // BURNIN/PROMETHEUS -> (BURNIN+STRUCTURED)/PROMETHEUS (30 -> 50)
  for (era = (int) BURNIN / PROMETHEUS; era < (int) (BURNIN + STRUCTURED) / PROMETHEUS; era++) {
    // GENEALOGY (with recombination and taking into account that duplicated chr have duplicated ancestor)
    cout << "running phaseII() era = " << era << '\n';
    genealogy(rho * BLOCKLENGTH, 1, (2 * k * BLOCKLENGTH));
    int prom=-1;
    prev = 0;
    pres = 1;
    bool skip = false;
    for (std::vector<fertility_info>::iterator it=fertility_list.begin(); it!=fertility_list.end(); ++it){
      int i = (*it).x;
      int t = (*it).y;
      if(prom!=t){
	if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
      }
      prom = t;
      if(skip == false){
	parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,i,t);// PARENT PICKING (with recombination)
	mutation(mu, i,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	if(IGCmatrix[i][t]==true && (i%2 == 0 )){
	  // WHC: if IGC happend, need to pick and mutate the partner chrom
	  // WHC: and this should only be done once, thus skip = true after
	  skip = true;
	  int otheri = i+1;
	  parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,otheri,t);// PARENT PICKING (with recombination)
	  mutation(mu, otheri,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	}
	conversion(kappa, t, i, pres, donorRatio, sameDifIGC);							
      }else {skip = false;}
    }
    // CALCULATE THE STATISTICS
    statistics(pres, does_print);
  }
}

void phaseIII(float k){
  bool does_print = false;
  bool run = false;
  for (era = (int) (BURNIN + STRUCTURED) / PROMETHEUS; era < (int) TIMELENGTH / PROMETHEUS; era++) {
    // GENEALOGY (all chr have the duplication, there are no two populations to take into account)
    cout << "running phaseIII() era = " << era << '\n';
    genealogy(rho * BLOCKLENGTH, 0, (2 * k * BLOCKLENGTH));
    int prom=-1;
    int prev = 0;
    int pres = 1;
    bool skip = false;
    for (std::vector<fertility_info>::iterator it=fertility_list.begin(); it!=fertility_list.end(); ++it){
      int i = (*it).x;
      int t = (*it).y;
      if(prom!=t){
	if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
      }
      prom = t;
      if(skip == false){
	parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,i,t);// PARENT PICKING (with recombination)
	mutation(mu, i,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	if(IGCmatrix[i][t]==true && (i%2 == 0)){
	  skip = true;
	  int otheri = i+1;
	  parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,otheri,t);// PARENT PICKING (with recombination)
	  mutation(mu, otheri,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	}
	conversion(kappa, t, i, pres, donorRatio, sameDifIGC);							
      }else {skip = false;}
    }
    // CALCULATE THE STATISTICS
    if(era == ((int) (TIMELENGTH/PROMETHEUS)-1)){
      does_print=true;
    }

    /* WHC: testing whether a mutational position has occured more than once for each block
    if (run == false) {
      for (int r = 0; r < MutCount; ++r) {
	cout << muttable[r].position << " " << muttable[r].block << endl;
      }
      run = true;
    }
    */
    statistics(pres, does_print);

  }
}

/*=====================================================================================================================*/
// WHC: my phaseIV() is a similar function to phaseII(), now

void phaseIV(int timeToFixation,int prev, int pres, float k){
  bool does_print = false;
  bool run = false;
  // Generating Fixation Trajectory in a diploid population
  int endTime = GenerateFixationTrajectory(STRUCTURED_2 + 1, timeToFixation);
  // Picking the first chromosome with the duplication (eva)
  int eva = (int) (rand() % (2 * N));
  for (int i = 0; i < 2 * N; i++) {
    duplicontent_2[0][i] = i;
    duplicontent_2[1][i] = i;
  }
  duplicontent_2[1][0] = eva;
  duplicontent_2[1][eva] = 0;

  // ================
  // WHC: up until here, duplication() need to be changed for phaseIV(); maybe consider writing a new duplication_2()?
  // Duplicate eva (create a new block with old mutations...)
  duplication_2(eva,pres,dupType); // WHC: the pres (present) here, will be used as prev (previous) in duplication(); because this is the start of PhaseII
  duFreq_2 = true;

  // BURNIN/PROMETHEUS -> (BURNIN+STRUCTURED)/PROMETHEUS (30 -> 50)
  // ==========================================================================================================================
  // WHC: stoped here
  
  //  for (era = (int) BURNIN / PROMETHEUS; era < (int) (BURNIN + STRUCTURED) / PROMETHEUS; era++) {
  for (era = (int) TIMELENGTH / PROMETHEUS; era < (int) (TIMELENGTH + STRUCTURED_2) / PROMETHEUS; era++) {
    // WHC: cout << "TIMELENGTH = " << TIMELENGTH << " " << "TIMELENGTH + STRUCTURED_2 = " << TIMELENGTH + STRUCTURED_2 << '\n';
    // GENEALOGY (with recombination and taking into account that duplicated chr have duplicated ancestor)
    //    genealogy(rho * BLOCKLENGTH, 1, (2 * k * BLOCKLENGTH));
    cout << "running phaseIV() era = " << era << '\n';
    genealogy_2(rho * BLOCKLENGTH, 1, (3 * k * BLOCKLENGTH));
    int prom=-1;
    prev = 0;
    pres = 1;
    bool skip = false;

    //    how_many_new_mutations_in_dup = 0;
    
    // WHC: I think it is good until here
    for (std::vector<fertility_info>::iterator it=fertility_list.begin(); it!=fertility_list.end(); ++it){
      // as long as PROMETHEUS is an even number, always end up being prov = 0, pres = 1
      int i = (*it).x;
      int t = (*it).y;
      if(prom!=t){
	if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
      }
      prom = t;
      if(skip == false){
	// WHC: I think parentpicking() works for phaseIV(); same is true for mutation()
	parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,i,t);// PARENT PICKING (with recombination)
	mutation(mu, i,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	if(IGCmatrix[i][t]==true && (i%2 == 0 )){
	  // WHC: if IGC happend, need to pick and mutate the partner chrom
	  // WHC: and this should only be done once, thus skip = true after
	  skip = true;
	  int otheri = i+1;
	  parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,otheri,t);// PARENT PICKING (with recombination)
	  mutation(mu, otheri,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	}
	conversion(kappa, t, i, pres, donorRatio, sameDifIGC);							
      }else {skip = false;}
    }
    // CALCULATE THE STATISTICS

    //    cout << "Then, there are " << how_many_new_mutations_in_dup << " introduced.\n\n";
    
    statistics(pres, does_print);

    /* ================================================================ */
    /* WHC: just to see how many chroms are carrying dup_2 */
    //    int counting_dup_2 = 0;
    //    for (int temp = 0; temp < 2 * N; ++temp) {
    //      if (pointer[pres][temp]->b == 5) {
    //	++counting_dup_2;
    //      }
    //    }
    //    cout << "The number of chroms that are carrying dup_2 = " << counting_dup_2 << '\n';
    /* WHC: the end of counting
    /* ================================================================ */

    /* WHC: testing for whether a muttalbe[] is added twice
    if (run == false) {
      for (int r = 0; r < MutCount; ++r) {
	cout << muttable[r].position << " " << muttable[r].block << endl;
      }
      run = true;
    }
    */

    
  }
}

/*=====================================================================================================================*/



/*=====================================================================================================================*/
// WHC: my phaseV() is a similar function to phaseII(), now

void phaseV(float k){
  bool does_print = false;
  bool run = false;
  for (era = (int) (TIMELENGTH + STRUCTURED_2) / PROMETHEUS; era < (int) (TIMELENGTH + STRUCTURED_2 + PHASE_V_LENGTH) / PROMETHEUS; era++) {

    cout << "Running phaseV era = " << era << '\n';
    // GENEALOGY (all chr have the duplication, there are no two populations to take into account)
    genealogy_2(rho * BLOCKLENGTH, 0, (3 * k * BLOCKLENGTH)); // WHC: genealogy_2(x, 0, x), the 0 means it is not in phaseIV()
    int prom=-1;
    int prev = 0;
    int pres = 1;
    bool skip = false;

    //    how_many_new_mutations_in_dup = 0;
    //    cout << "The length of this fertility_list is " << fertility_list.size() << endl;
    for (std::vector<fertility_info>::iterator it=fertility_list.begin(); it!=fertility_list.end(); ++it){
      int i = (*it).x;
      int t = (*it).y;
      if(prom!=t){
	if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
      }
      prom = t;
      if(skip == false){
	parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,i,t);// PARENT PICKING (with recombination)

	// WHC: use phaseVI's function to test
	
	//	parentpicking_for_phaseVI(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,i,t);
	//	mutation(mu, i,pres);// MUTATION and CONVERSION (for each fertile chromosome)

	//mutation_for_phaseVI(mu, i,pres);
	mutation(mu, i,pres);
	
	if(IGCmatrix[i][t]==true && (i%2 == 0)){
	  skip = true;
	  int otheri = i+1;
	  parentpicking(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,otheri,t);// PARENT PICKING (with recombination)
	  mutation(mu, otheri,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	  //mutation_for_phaseVI(mu, otheri,pres);
	}
	conversion(kappa, t, i, pres, donorRatio, sameDifIGC);							
      }else {skip = false;}
    }
    // CALCULATE THE STATISTICS
    if(era == ((int) (TIMELENGTH/PROMETHEUS)-1)){
      does_print=true;
    }

    //    cout << "Then, there are " << how_many_new_mutations_in_dup << " introduced.\n\n";
    
    /* WHC: testing, as before, to see if muttable[] has redundency
    if (run == false) {
      for (int r = 0; r < MutCount; ++r) {
	cout << muttable[r].position << " " << muttable[r].block << endl;
      }
      run = true;
    }
    */
    statistics(pres, does_print);
    //    statistics_for_phaseVI(pres, does_print);
    
  }
}

/*=====================================================================================================================*/

/*=====================================================================================================================*/
// WHC: my phaseVI() is a similar function to phaseIV(), but uses its own Generate_phaseVI_trajectory() function

void phaseVI(int prev, int pres, float k) {
  bool does_print = false;
  bool run = false;
  // Generating Fixation Trajectory in a diploid population
  Generate_phaseVI_trajectory(2); // WHC: trajectory (absolute number of chroms carrying 4 blocks instead of 5)

  // WHC: the Generate_phaseVI_trajectory() will generate #chroms that lose dup_1
  // WHC: however, as crossover will change it...
  // WHC: so, lose_of_duplicontent[previous][] must be re-sorted

  
  // Picking the first chromosome with the lose of dup_1 (eva)
  int eva = (int) (rand() % (2 * N));
  for (int i = 0; i < 2 * N; i++) {
    lose_of_duplicontent[0][i] = i;
    lose_of_duplicontent[1][i] = i;
      //    duplicontent_2[0][i] = i;
      //    duplicontent_2[1][i] = i;
  }
  //  duplicontent_2[1][0] = eva;
  //  duplicontent_2[1][eva] = 0;

  lose_of_duplicontent[1][0] = eva;
  lose_of_duplicontent[1][eva] = 0;
  
  // ================
  // WHC: lose_of_duplication() is for losing one dup_1 block in one chrom (i.e. eva)
  // Duplicate eva (create a new block with old mutations...)
  //  duplication_2(eva,pres,dupType); // WHC: the pres (present) here, will be used as prev (previous) in duplication(); because this is the start of PhaseII



  // WHC: to test what's wrong

    lose_of_duplication(eva, pres); // WHC: as will lose dup_1 block from eva chrom, there is no need for the dupType parameter




  
  loseFreq = true;

  // BURNIN/PROMETHEUS -> (BURNIN+STRUCTURED)/PROMETHEUS (30 -> 50)
  // ==========================================================================================================================
  // WHC: stoped here
  
  //  for (era = (int) TIMELENGTH / PROMETHEUS; era < (int) (TIMELENGTH + STRUCTURED_2) / PROMETHEUS; era++) {
  for (era = (int) (TIMELENGTH + STRUCTURED_2 + PHASE_V_LENGTH) / PROMETHEUS; era < (int) (TIMELENGTH + STRUCTURED_2 + PHASE_V_LENGTH + PHASE_VI_LENGTH) / PROMETHEUS; era++) {

    cout << "Running phaseVI, era = " << era << '\n';
    
    // WHC: cout << "TIMELENGTH = " << TIMELENGTH << " " << "TIMELENGTH + STRUCTURED_2 = " << TIMELENGTH + STRUCTURED_2 << '\n';
    // GENEALOGY (with recombination and taking into account that duplicated chr have duplicated ancestor)
    //    genealogy(rho * BLOCKLENGTH, 1, (2 * k * BLOCKLENGTH));


    genealogy_2(rho * BLOCKLENGTH, 2, (3 * k * BLOCKLENGTH)); // WHC: 2 means it is in phaseVI(), not phaseIV()


    //    genealogy_2(rho * BLOCKLENGTH, 0, (3 * k * BLOCKLENGTH));
    
    int prom=-1;
    prev = 0;
    pres = 1;
    bool skip = false;

    //    how_many_new_mutations_in_dup = 0;
    
    // WHC: I think it is good until here
    for (std::vector<fertility_info>::iterator it=fertility_list.begin(); it!=fertility_list.end(); ++it){
      // as long as PROMETHEUS is an even number, always end up being prov = 0, pres = 1
      int i = (*it).x;
      int t = (*it).y;
      if(prom!=t){
	if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
      }
      prom = t;
      if(skip == false){
	// WHC: I think parentpicking() works for phaseIV(); same is true for mutation()

	// WHC: has to decide which blcok does the crossover happen, before decidiing #blocks for child chrom!!!
	
	parentpicking_for_phaseVI(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,i,t);// PARENT PICKING (with recombination)
	//	mutation(mu, i,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	mutation_for_phaseVI(mu, i, pres);
	if(IGCmatrix[i][t]==true && (i%2 == 0 )){
	  // WHC: if IGC happend, need to pick and mutate the partner chrom
	  // WHC: and this should only be done once, thus skip = true after
	  skip = true;
	  int otheri = i+1;
	  parentpicking_for_phaseVI(crossoverBegin, crossoverEnd, crossoverRatio, numHS,prev,pres,otheri,t);// PARENT PICKING (with recombination)
	  mutation_for_phaseVI(mu, otheri,pres);// MUTATION and CONVERSION (for each fertile chromosome)
	}
	conversion_for_phaseVI(kappa, t, i, pres, donorRatio, sameDifIGC);							
      } else {skip = false;}
    }
    // CALCULATE THE STATISTICS

    // WHC: testing for muttable[] redudency
    /*
    if (run == false) {
      for (int r = 0; r < MutCount; ++r) {
	if (muttable[r].block == 1) {
	  //	cout << muttable[r].position << " " << muttable[r].block << endl;
	  cout << muttable[r].position << endl;
	}
      }
      run = true;
    }
    */

    //    cout << "Then, there are " << how_many_new_mutations_in_dup << " introduced.\n\n";    
    statistics_for_phaseVI(pres, does_print);

    
    // WHC: for printing block_1's mutation positions
    /*
    for (int i = 0; i < 2 * N; ++i) {
      cout << "this is for individual " << i << '\n';
      for (int k = 0; k < pointer[pres][i]->mpb[1]; ++k) {
	cout << k << '\t' << pointer[pres][i]->mutation[1][k] << '\n';
      }
    }
    */

    
    /* ================================================================ */
    /* WHC: just to see how many chroms are carrying dup_2 */
    /*
    int counting_losing_dup_1 = 0, counting_having_dup_1 = 0;
        for (int temp = 0; temp < 2 * N; ++temp) {
          if (pointer[pres][temp]->b == 4) {
	    ++counting_losing_dup_1;
          } else if (pointer[pres][temp]->b == 5) {
	    ++counting_having_dup_1;
	  } else { cout << "WRONG HERE!\n"; exit(0); }
        }
        cout << "The number of chroms that are not carrying dup_1 = " << counting_losing_dup_1 << '\n';
	cout << "The number of chroms that are carrying dup_1 = " << counting_having_dup_1 << '\n';
    */
    /* WHC: the end of counting
    /* ================================================================ */
    
  }
}


/*=====================================================================================================================*/





//////////////////////////////
////  OPEN & CLOSE FILES  ////
//////////////////////////////

void open_files() { // Opens write-on files
  int j, o;
  stringstream ss;

  //ss << "auxx_" << letter << ".dat" << endl;
  //ss >> str;
  //auxx.open(str.c_str());
  //ss.str("");

  ss << "profile_" << letter << ".dat" << endl;
  ss >> str;
  profile.open(str.c_str());
  ss.str("");

  for (j = 0; j < B; j++) {
    ss << "mutations_new_" << j << "_" << letter << ".dat" << endl;
    ss >> str;
    mutationsNewFile[j].open(str.c_str());
    ss.str("");
  }

  for (j = 0; j < B; j++) {
    ss << "SFS_" << j << "_" << letter << ".dat" << endl;
    ss >> str;
    SFS[j].open(str.c_str());
    ss.str("");
  }

  // WHC: create files for inter-block pairwise divergence

  interblock_divergence[0].open("interblock_divergence_between_[0][3].dat"); // WHC: between original and dup_1
  interblock_divergence[1].open("interblock_divergence_between_[0][5].dat"); // WHC: between original and dup_2
  interblock_divergence[2].open("interblock_divergence_between_[3][5].dat"); // WHC: between dup_1 and dup_2

  for (o = 0; o < numofsamples; o++) {
    for (j = 0; j < B; j++) {
      ss << "samplepi" << j << "[" << sampleN[o] << "]_" << letter << ".dat" << endl;
      ss >> str;
      samplefile[j][o]->open(str.c_str());
      ss.str("");

      ss << "sampleS" << j << "[" << sampleN[o] << "]_" << letter << ".dat" << endl;
      ss >> str;
      samplefile[j][o][1].open(str.c_str());
      ss.str("");
    }
  }
}

void close_files() { // Closes write-on files
  int j, o;

  profile.close();
  //auxx.close();
  for (j = 0; j < B; j++) { mutationsNewFile[j].close();}
  for (j = 0; j < B; j++) {
    SFS[j].close();
    for (o = 0; o < numofsamples; o++) {
      samplefile[j][o]->close();
      samplefile[j][o][1].close();
    }
  }

  // WHC: close files for inter-block divergence
  for (int k = 0; k < 3; ++k) {
    interblock_divergence[k].close();
  }
}

////////////////////////////////////////////
////  GENEALOGY & STRUCTURED GENEALOGY  ////
////////////////////////////////////////////

void genealogy(float probability, int strornot, float IGCprobability) { // Generates the genealogy based on the FixationTrajectory (with RECOMBINATION)
  // Determine recombination processes (with probability "probability")
  for (int tt=0 ; tt < PROMETHEUS ; tt++){
    for (int i=0 ; i < 2*N ; i++){
      recombimatrix[i][tt] = false;
      IGCmatrix[i][tt] = false;
    }
  }

  // STRUCTURED GENEALOGY (the population is subdivided in 2: one that carries the duplication, built according to trajectime, and the other not carrying the duplication
  if (strornot == 1){
    //cout << "entra a strornot\n";
    // Pick a parent from the corresponding duplicated/non-duplicated population
    int present=0;
    int previous=1;
    for (int tt=0 ; tt < PROMETHEUS ; tt++){
      // trajectime is the time in which we are in fixationTrajectory[] array (used to know the number of chr that carry the dup at each time)
      int trajectime = PROMETHEUS*(era-((int)BURNIN/PROMETHEUS))+tt+1;
      // Randomly mix all the chromosomes of the present generation
      // WHC: as the next step (choosing a father with dup ramdomly always pick the first several chroms, this step is for randomly picking chroms for having dup
      
      for(int i=0 ; i < 2*N ; i++){
	int val = (int) (rand()%(2*N));
	int temp = duplicontent[present][i];
	// Duplicontent is an array that indicates which chr carry the duplication (with present and previous lines)
	duplicontent[present][i] = duplicontent[present][val];
	duplicontent[present][val] = temp;
      }
      // For all the chr that have to have the duplication at this time...
      for(int i=0 ; i < fixationTrajectory[trajectime] ; i++){
	// Choose a father randomly (from the previous duplicated population)
	// WHC: as present chroms which will have dup have been chosen before, this is for randomly picking a father
	// WHC: remember, i = 0, 1, 2, ...
	int val=(int) (rand()%(fixationTrajectory[trajectime-1]));
	ancestry[duplicontent[present][i]][tt] = duplicontent[previous][val];
      }
      // For all the chr that have not to have the duplication at this time...
      for(int i=fixationTrajectory[trajectime] ; i < 2*N ; i++){
	// Choose a father randomly (from the previous non-duplicated population)
	int val = (int) (rand()%(2*N-fixationTrajectory[trajectime-1])) + fixationTrajectory[trajectime-1];
	ancestry[duplicontent[present][i]][tt] = duplicontent[previous][val];
      }
      if (previous==1) {previous=0;  present=1;} else {previous=1;  present=0;}
    }
  }
  // NORMAL GENEALOGY (no 2 populations Duplicated/NoDuplicated)
  else if (strornot == 0){
    for (int tt=0 ; tt < PROMETHEUS ; tt++){
      for (int i=0 ; i < 2*N ; i++){
	int val=(int) (rand()%(2*N));
	ancestry[i][tt] = val; //   WITHOUT presLOSS OF GENERALITY AN INDIV. CAN MATE WITH HIMSELF
      }
    }
  }
  // TRACING BACK THE GENEALOGY AND BUILDING FERTILITY MATRIX
  fertility = fertility_ini;
  fertility_list.clear();
  for (int i=0 ; i < 2*N ; i++){
    fertility_info fi;
    fi.x = i;
    fi.y = (PROMETHEUS-1);
    fertility_list.push_back(fi);
  }
  for (int tt=PROMETHEUS-1 ; tt > 0 ; tt--){
    for (int i=0 ; i < 2*N ; i++){
      if (fertility[i][tt] == true){
	if(fertility[ancestry[i][tt]][tt-1]==false){
	  fertility_info fi;
	  fi.x=(ancestry[i][tt]);
	  fi.y=(tt-1);
	  fertility_list.push_back(fi);
	}
	fertility[ancestry[i][tt]][tt-1] = true;
	float p = rand() / ((float) RAND_MAX + 1);// Determine recombination processes (with probability "probability")
	if (p < probability){
	  recombimatrix[i][tt] = true;
	  if(ancestry[i][tt]%2==0) {
	    if(fertility[ancestry[i][tt]+1][tt-1]==false){
	      fertility_info fi;
	      fi.x=(ancestry[i][tt]+1);
	      fi.y=(tt-1);
	      fertility_list.push_back(fi);
	    }
	    fertility[ancestry[i][tt]+1][tt-1]=true;
	  }
	  else {
	    if(fertility[ancestry[i][tt]-1][tt-1]==false){
	      fertility_info fi;
	      fi.x=(ancestry[i][tt]-1);
	      fi.y=(tt-1);
	      fertility_list.push_back(fi);
	    }
	    fertility[ancestry[i][tt]-1][tt-1]=true;
	  }
	}
			
	if((sameDifIGC!=1)){
	  p = rand() / ((float) RAND_MAX + 1);// Determine IGC processes (with probability "IGCprobability"); WHC: probabily error, should be "Determin IGC processes"; also, it is determining IGC process between 2 chroms, instead of on itself

	  // WHC: the (1 - sameDifIGC), is due to the fact there is another p < (2 * probability * BLOCKLENGTH * sameDifIGC) in void conversion(); this one here is to only test whether IGC comes from different chromosomes, if so, need to record its genealogy; if not, then it either doesn't have IGC, or it has IGC within the same chromosome, which doesn't require genealogy for its partner chromosome
	  
	  if (p < (IGCprobability * (1-sameDifIGC))){
	    IGCmatrix[i][tt] = true;
	    if(i%2==0) {
	      // WHC: here, fertility[i+1][tt] will be true, which means [i+1] will be tested for IGC & recombination
	      // WHC: however, if i%2 == 1, which means its partner [i-1] may have fertility[i-1][tt] == false, which means
	      // WHC: it is NOT tested for recombination
	      // WHC: that is why the else {} tests recombination for it
	      if(fertility[i+1][tt]==false){
		fertility_info fi;
		fi.x=(i+1);
		fi.y=(tt);
		fertility_list.push_back(fi);
	      }
	      fertility[i+1][tt]=true;
	    }else{
	      // WHC: same algorithm as in void phaseII() -- scan from 1st chrom to 2*N th chrom; if i%2 == 0, then record (i+1)'s (its partner that hasn't been scanned yet) genealogy; if i%2 != 0, which means its previous partner is already scanned
	      
	      if(fertility[i-1][tt]==false){
		fertility_info fi;
		fi.x=(i-1);
		fi.y=(tt);
		fertility_list.push_back(fi);
							
		fertility[i-1][tt]=true;
								
		if(fertility[ancestry[i-1][tt]][tt-1]==false){
		  fertility_info fi;
		  fi.x=(ancestry[i-1][tt]);
		  fi.y=(tt-1);
		  fertility_list.push_back(fi);
		}
		fertility[ancestry[i-1][tt]][tt-1] = true;
							 
		p = rand() / ((float) RAND_MAX + 1);// Determine recombination processes (with probability "probability")
		if (p < probability){
		  recombimatrix[i-1][tt] = true;
		  if(ancestry[i-1][tt]%2==0) {
		    if(fertility[ancestry[i-1][tt]+1][tt-1]==false){
		      fertility_info fi;
		      fi.x=(ancestry[i-1][tt]+1);
		      fi.y=(tt-1);
		      fertility_list.push_back(fi);
		    }
		    fertility[ancestry[i-1][tt]+1][tt-1]=true;
		  }else {
		    if(fertility[ancestry[i-1][tt]-1][tt-1]==false){
		      fertility_info fi;
		      fi.x=(ancestry[i-1][tt]-1);
		      fi.y=(tt-1);
		      fertility_list.push_back(fi);
		    }
		    fertility[ancestry[i-1][tt]-1][tt-1]=true;
		  }
		}
		
		// WHC: this p = rand() here, is due to, if i % 2 == 0, then (i+1) will be added to fertility_list
		// WHC: which means that it will NOT be tested 
		
		p = rand() / ((float) RAND_MAX + 1);// Determine recombination processes (with probability "probability")
		if (p < (IGCprobability * (1-sameDifIGC))){
		  IGCmatrix[i-1][tt] = true;
		}
	      }
	    }
	  }
	}	
      }
    }
  }
	
  std::stable_sort (fertility_list.begin(), fertility_list.end(), sortx);
  std::stable_sort (fertility_list.begin(), fertility_list.end(), sorty);
		
}

/* ======================================================================================================================= */
// WHC: genealogy_2() for phaseIV()

void genealogy_2(float probability, int strornot, float IGCprobability) { // Generates the genealogy based on the FixationTrajectory (with RECOMBINATION)
  // Determine recombination processes (with probability "probability")
  for (int tt=0 ; tt < PROMETHEUS ; tt++){
    for (int i=0 ; i < 2*N ; i++){
      recombimatrix[i][tt] = false;
      IGCmatrix[i][tt] = false;
    }
  }

  // STRUCTURED GENEALOGY (the population is subdivided in 2: one that carries the duplication, built according to trajectime, and the other not carrying the duplication
  
  // WHC: this strornot == 1 means in phaseIV()
  if (strornot == 1){
    //cout << "entra a strornot\n";
    // Pick a parent from the corresponding duplicated/non-duplicated population
    int present=0;
    
    // WHC: previous = 1 is last present = 1; now becomes previous generation
    int previous=1;
    for (int tt=0 ; tt < PROMETHEUS ; tt++){
      // trajectime is the time in which we are in fixationTrajectory[] array (used to know the number of chr that carry the dup at each time)
      int trajectime = PROMETHEUS*(era-((int)TIMELENGTH/PROMETHEUS))+tt+1;
      // Randomly mix all the chromosomes of the present generation
      // WHC: as the next step (choosing a father with dup ramdomly always pick the first several chroms, this step is for randomly picking chroms for having dup
      
      for(int i=0 ; i < 2*N ; i++){
	int val = (int) (rand()%(2*N));
	int temp = duplicontent_2[present][i];
	// Duplicontent is an array that indicates which chr carry the duplication (with present and previous lines)
	duplicontent_2[present][i] = duplicontent_2[present][val];
	duplicontent_2[present][val] = temp;
      }
      // For all the chr that have to have the duplication at this time...

      // WHC: need to generate a fixationTrajectory_2[] for phaseIV() ?
      for(int i=0 ; i < fixationTrajectory[trajectime] ; i++){
	// Choose a father randomly (from the previous duplicated population)
	// WHC: as present chroms which will have dup have been chosen before, this is for randomly picking a father
	// WHC: remember, i = 0, 1, 2, ...
	int val=(int) (rand()%(fixationTrajectory[trajectime-1]));
	
	// WHC: FOUND a MISTAKE! should be duplicatent_2[]!!
	// ancestry[duplicontent[present][i]][tt] = duplicontent[previous][val];
	ancestry[duplicontent_2[present][i]][tt] = duplicontent_2[previous][val];
      }
      // For all the chr that have not to have the duplication at this time...
      for(int i=fixationTrajectory[trajectime] ; i < 2*N ; i++){
	// Choose a father randomly (from the previous non-duplicated population)
	int val = (int) (rand()%(2*N-fixationTrajectory[trajectime-1])) + fixationTrajectory[trajectime-1];
	ancestry[duplicontent_2[present][i]][tt] = duplicontent_2[previous][val];
      }
      if (previous==1) {previous=0;  present=1;} else {previous=1;  present=0;}
    }
  } else if (strornot == 0){  // NORMAL GENEALOGY (no 2 populations Duplicated/NoDuplicated); in phaseV()
    for (int tt=0 ; tt < PROMETHEUS ; tt++){
      for (int i=0 ; i < 2*N ; i++){
	int val=(int) (rand()%(2*N));
	ancestry[i][tt] = val; //   WITHOUT presLOSS OF GENERALITY AN INDIV. CAN MATE WITH HIMSELF
      }
    }
  } else if (strornot == 2) {	// WHC: means it is in phaseVI(), not in phaseIV(); strornot == 0 is not necessary in genealogy_2()
    /* this is for generating the genealogy for phaseVI, the losing of dup_1 for some chroms */

    int present=0;
    
    // WHC: previous = 1 is last present = 1; now becomes previous generation
    int previous=1;
    for (int tt=0 ; tt < PROMETHEUS ; tt++){
      // trajectime is the time in which we are in fixationTrajectory[] array (used to know the number of chr that carry the dup at each time)
      // WHC: as phaseVI_trajectory[] 
      int trajectime = PROMETHEUS*(era-((int) (TIMELENGTH + STRUCTURED_2 + PHASE_V_LENGTH) / PROMETHEUS))+tt+1;
      

      
      // Randomly mix all the chromosomes of the present generation
      // WHC: as the next step (choosing a father with dup ramdomly always pick the first several chroms, this step is for randomly picking chroms for having dup
      
      for(int i=0 ; i < 2*N ; i++){
	int val = (int) (rand()%(2*N));
	int temp = lose_of_duplicontent[present][i];
	// Duplicontent is an array that indicates which chr carry the duplication (with present and previous lines)
	lose_of_duplicontent[present][i] = lose_of_duplicontent[present][val];
	lose_of_duplicontent[present][val] = temp;
      }
      // For all the chr that have to have the duplication at this time...

      // WHC: need to generate a fixationTrajectory_2[] for phaseIV() ?
      for(int i=0 ; i < phaseVI_trajectory[trajectime] ; i++){
	// Choose a father randomly (from the previous duplicated population)
	// WHC: as present chroms which will have dup have been chosen before, this is for randomly picking a father
	// WHC: remember, i = 0, 1, 2, ...
	int val=(int) (rand()%(phaseVI_trajectory[trajectime-1]));

	// WHC: also, phaseVI_trajectory[0] == 1!!! as is the case for fixationTrajectory[0] = 1
	
	// ancestry[duplicontent[present][i]][tt] = duplicontent[previous][val];
	ancestry[lose_of_duplicontent[present][i]][tt] = lose_of_duplicontent[previous][val];
      }
      
      // For all the chr that have not to have the duplication at this time...
      for(int i=phaseVI_trajectory[trajectime] ; i < 2*N ; i++){
	// Choose a father randomly (from the previous non-duplicated population)
	int val = (int) (rand()%(2*N-phaseVI_trajectory[trajectime-1])) + phaseVI_trajectory[trajectime-1];
	ancestry[lose_of_duplicontent[present][i]][tt] = lose_of_duplicontent[previous][val];
      }
      if (previous==1) {previous=0;  present=1;} else {previous=1;  present=0;}
    }

    /* the end of this clause */
  }
  // TRACING BACK THE GENEALOGY AND BUILDING FERTILITY MATRIX
  fertility = fertility_ini;
  fertility_list.clear();
  for (int i=0 ; i < 2*N ; i++){
    fertility_info fi;
    fi.x = i;
    fi.y = (PROMETHEUS-1);
    fertility_list.push_back(fi);
  }
  for (int tt=PROMETHEUS-1 ; tt > 0 ; tt--){
    for (int i=0 ; i < 2*N ; i++){
      if (fertility[i][tt] == true){
	if(fertility[ancestry[i][tt]][tt-1]==false){
	  fertility_info fi;
	  fi.x=(ancestry[i][tt]);
	  fi.y=(tt-1);
	  fertility_list.push_back(fi);
	}
	fertility[ancestry[i][tt]][tt-1] = true;
	float p = rand() / ((float) RAND_MAX + 1);// Determine recombination processes (with probability "probability")
	if (p < probability){
	  recombimatrix[i][tt] = true;
	  if(ancestry[i][tt]%2==0) {
	    if(fertility[ancestry[i][tt]+1][tt-1]==false){
	      fertility_info fi;
	      fi.x=(ancestry[i][tt]+1);
	      fi.y=(tt-1);
	      fertility_list.push_back(fi);
	    }
	    fertility[ancestry[i][tt]+1][tt-1]=true;
	  }
	  else {
	    if(fertility[ancestry[i][tt]-1][tt-1]==false){
	      fertility_info fi;
	      fi.x=(ancestry[i][tt]-1);
	      fi.y=(tt-1);
	      fertility_list.push_back(fi);
	    }
	    fertility[ancestry[i][tt]-1][tt-1]=true;
	  }
	}
			
	if((sameDifIGC!=1)){
	  p = rand() / ((float) RAND_MAX + 1);// Determine IGC processes (with probability "IGCprobability"); WHC: probabily error, should be "Determine IGC processes"; also, it is determining IGC process between 2 chroms, instead of on itself

	  // WHC: the (1 - sameDifIGC), is due to the fact there is another p < (2 * probability * BLOCKLENGTH * sameDifIGC) in void conversion(); this one here is to only test whether IGC comes from different chromosomes, if so, need to record its genealogy; if not, then it either doesn't have IGC, or it has IGC within the same chromosome, which doesn't require genealogy for its partner chromosome
	  
	  if (p < (IGCprobability * (1-sameDifIGC))){
	    IGCmatrix[i][tt] = true;
	    if(i%2==0) {
	      // WHC: here, fertility[i+1][tt] will be true, which means [i+1] will be tested for IGC & recombination
	      // WHC: however, if i%2 == 1, which means its partner [i-1] may have fertility[i-1][tt] == false, which means
	      // WHC: it is NOT tested for recombination
	      // WHC: that is why the else {} tests recombination for it
	      if(fertility[i+1][tt]==false){
		fertility_info fi;
		fi.x=(i+1);
		fi.y=(tt);
		fertility_list.push_back(fi);
	      }
	      fertility[i+1][tt]=true;
	    }else{
	      // WHC: same algorithm as in void phaseII() -- scan from 1st chrom to 2*N th chrom; if i%2 == 0, then record (i+1)'s (its partner that hasn't been scanned yet) genealogy; if i%2 != 0, which means its previous partner is already scanned
	      
	      if(fertility[i-1][tt]==false){
		fertility_info fi;
		fi.x=(i-1);
		fi.y=(tt);
		fertility_list.push_back(fi);
							
		fertility[i-1][tt]=true;
								
		if(fertility[ancestry[i-1][tt]][tt-1]==false){
		  fertility_info fi;
		  fi.x=(ancestry[i-1][tt]);
		  fi.y=(tt-1);
		  fertility_list.push_back(fi);
		}
		fertility[ancestry[i-1][tt]][tt-1] = true;
							 
		p = rand() / ((float) RAND_MAX + 1);// Determine recombination processes (with probability "probability")
		if (p < probability){
		  recombimatrix[i-1][tt] = true;
		  if(ancestry[i-1][tt]%2==0) {
		    if(fertility[ancestry[i-1][tt]+1][tt-1]==false){
		      fertility_info fi;
		      fi.x=(ancestry[i-1][tt]+1);
		      fi.y=(tt-1);
		      fertility_list.push_back(fi);
		    }
		    fertility[ancestry[i-1][tt]+1][tt-1]=true;
		  }else {
		    if(fertility[ancestry[i-1][tt]-1][tt-1]==false){
		      fertility_info fi;
		      fi.x=(ancestry[i-1][tt]-1);
		      fi.y=(tt-1);
		      fertility_list.push_back(fi);
		    }
		    fertility[ancestry[i-1][tt]-1][tt-1]=true;
		  }
		}
		
		// WHC: this p = rand() here, is due to, if i % 2 == 0, then (i+1) will be added to fertility_list
		// WHC: which means that it will NOT be tested 
		
		p = rand() / ((float) RAND_MAX + 1);// Determine recombination processes (with probability "probability")
		if (p < (IGCprobability * (1-sameDifIGC))){
		  IGCmatrix[i-1][tt] = true;
		}
	      }
	    }
	  }
	}	
      }
    }
  }
	
  std::stable_sort (fertility_list.begin(), fertility_list.end(), sortx);
  std::stable_sort (fertility_list.begin(), fertility_list.end(), sorty);
		
}

// WHC: genealogy_2() for phaseIV()
/* ======================================================================================================================= */

////////////////////////////////////////////////////////////
////  PARENTPICKING, DUPLICATION, MUTATION, CONVERSION  ////
////////////////////////////////////////////////////////////

void parentpicking(int crossBegin[maxNumOfHS], int crossEnd[maxNumOfHS], float crossRatio[maxNumOfHS], int numCrossRegions, int prev, int pres, int i, int t) {
  int j, k, junctionBlock, defHS, HS;
  int father, partner, junction, vald0, valr0, childblocks, minblock, recTract, end, beg;
  bool success;
  chrom * chr;
  float p;
  double num;

  // PARENT-PICKING
  father = ancestry[i][t];
  // RECOMBINATION
  if (recombimatrix[i][t] == true) {
    childblocks = pointer[prev][father]->b;
    if (ancestry[i][t] % 2 == 0) {
      partner = ancestry[i][t] + 1;
    } else {
      partner = ancestry[i][t] - 1;
    }
    chr = pointer[pres][i];
    chr->b = childblocks;
    minblock = minim(pointer[prev][father]->b, pointer[prev][partner]->b);
    //DETERMINE THE HOTSPOT WHERE CROSSOVER WILL OCCUR

    // WHC: Should always consider crossover after deletion happened
    // WHC: probabily will just write a parentpicking() specifically for that phase...
    // WHC: otherwise I think this is workable for now; when a block that doesn't exist happen, just no crossover and copy parental chrom
    p = rand() / ((float) RAND_MAX + 1);
    defHS = numCrossRegions-1;
    num = 1;
    HS = numCrossRegions-1;
    success = 0;
    while(HS >= 0 && success==0){
      num = num-crossRatio[HS];
      if (p < num){
	defHS = HS-1;
      }
      else {
	success=1;
      }
      HS--;
    }
    end = crossEnd[defHS];
    beg = crossBegin[defHS];
    if (((minblock * BLOCKLENGTH) < end) && ((minblock * BLOCKLENGTH) > beg)) {
      end = minblock * BLOCKLENGTH;
    }
    if (((minblock * BLOCKLENGTH) >= end) && ((minblock * BLOCKLENGTH) > beg)) {
      // If neither father and partner have mutations in any block
      /*      if ((pointer[prev][father]->mpb[0] == 0) && (pointer[prev][father]->mpb[1] == 0) &&
	  (pointer[prev][father]->mpb[2] == 0) && (pointer[prev][partner]->mpb[0] == 0) &&
	  (pointer[prev][partner]->mpb[1] == 0) && (pointer[prev][partner]->mpb[2] == 0)) {
      */
      if ((pointer[prev][father]->mpb[0] == 0) && (pointer[prev][father]->mpb[1] == 0) &&
	  (pointer[prev][father]->mpb[2] == 0) && (pointer[prev][father]->mpb[3] == 0) &&
	  (pointer[prev][father]->mpb[4] == 0) && (pointer[prev][father]->mpb[5] == 0) &&
	  (pointer[prev][partner]->mpb[0] == 0) && (pointer[prev][partner]->mpb[1] == 0) &&
	  (pointer[prev][partner]->mpb[2] == 0) && (pointer[prev][partner]->mpb[3] == 0) &&
	  (pointer[prev][partner]->mpb[4] == 0) && (pointer[prev][partner]->mpb[5] == 0)) {
      for (j = 0; j < chr->b; j++) {
	chr->mpb[j] = 0;
	}
      }// If at least one of them have one or more mutations in one or more blocks
      else {

	recTract = end - beg;
	junction = (int) (rand() % (recTract));
	junction += beg;
	junctionBlock = (int) (junction / BLOCKLENGTH); // block where junction fell
	//	   cout << junctionBlock << " ";
	// WHC: is the number of mutations before junction point
	valr0 = location(junction - junctionBlock*BLOCKLENGTH, prev, partner, junctionBlock); // junction location in the block "junctionBlock" of partner
	vald0 = location(junction - junctionBlock*BLOCKLENGTH, prev, father, junctionBlock); // junction location in the block "junctionBlock" of father

	//COPY THE INFO FROM PARTNER (BLOCK 0)
	for (j = 0; j < junctionBlock; j++) {
	  chr->mpb[j] = pointer[prev][partner]->mpb[j];
	  for (k = 0; k < pointer[prev][partner]->mpb[j]; k++) {
	    chr->mutation[j][k] = pointer[prev][partner]->mutation[j][k];
	  }
	}

	//COPY INFO IN JUNCTION BLOCK (BLOCK 1)
	//FROM PARTNER

	// WHC: copy mutation from [partner] until junction point
	for (k = 0; k < valr0; k++) {
	  chr->mutation[junctionBlock][k] = pointer[prev][partner]->mutation[junctionBlock][k];
	}
	//FROM FATHER
	
	// WHC: copy mutation from [father] from junction point on
	for (k = 0; k < (pointer[prev][father]->mpb[junctionBlock] - vald0); k++) {
	  chr->mutation[junctionBlock][k + valr0] = pointer[prev][father]->mutation[junctionBlock][k + vald0];
	}
	//ESTABLISH MPB
	chr->mpb[junctionBlock] = valr0 + (pointer[prev][father]->mpb[junctionBlock] - vald0);

	//COPY INFO FROM BLOCK 2 IF PRESENT

	// WHC: WRONG??? forgot copying some information???

	
	// WHC: as from junction point on, all mutation info comes from father...
	// WHC: (also, previous step only copy mutation info with the junctionBlock!!!
	for (j = junctionBlock + 1; j < chr->b; j++) {
	  if (j < pointer[prev][father]->b) { //j counts 0,1,2 but .b counts 1,2,3
	    chr->mpb[j] = pointer[prev][father]->mpb[j];
	    for (k = 0 ; k < pointer[prev][father]->mpb[j] ; k++) {
	      chr->mutation[j][k] = pointer[prev][father]->mutation[j][k];
	    }
	  } else {
	    chr->mpb[j] = 0;
	  }
	}

	// WHC: IMPORTANT!!! Set mpb[4] == 0 if chr->b == 3
	//	if (chr->b == 3) {
	if (chr->b == 4) {
	  chr->mpb[5] = 0;
	  chr->mpb[4] = 0;
	  //	} else if (chr->b == 2) {
	} else if (chr->b == 2) {
	  chr->mpb[2] = 0;
	  chr->mpb[3] = 0;
	  chr->mpb[4] = 0;
	  chr->mpb[5] = 0;
	}
      }
    } else if (((minblock * BLOCKLENGTH) < end) && ((minblock * BLOCKLENGTH) <= beg)) {
      // WHC: this steps looks like a "lazy" step; (NO! I was wrong! Imagine if a hot block on block 3 with 100% proposion!
      // WHC: (then, if we don't have block 3, we don't have recombination!)
      copychr(prev, father, pres, i);
    }
  }// NO RECOMBINATION
  else {
    copychr(prev, father, pres, i);
  }

}

/* ================================================================ */
// WHC: parentpicking_for_phaseVI() is to deal with when some chroms have 4 blocks, some have 5
// WHC: which needs to pick crossover positon (and partner) before deciding child chrom's #blocks

void parentpicking_for_phaseVI(int crossBegin[maxNumOfHS], int crossEnd[maxNumOfHS], float crossRatio[maxNumOfHS], int numCrossRegions, int prev, int pres, int i, int t) {
  // WHC: only consider the situation where a chrom either has 4 or 5 blocks; and crossover can only happen on common blocks

  /************************************************************************************************************************
    WHC: IMPORTANT: for crossover between #block = 5 && # blocks = 4 only: (chr->b == #father's blocks
    if corssover happens in HS0 (block0, block1), then to make chr->b = 4, before junction pointer, all information comes from partner;
    else if crossover happens in HS1 or HS2, then before junction pointer, all informations comes from father.
    this is just to make sure, information about dup_1 (have or lose) comes from father
  ************************************************************************************************************************/
  
  int j, k, junctionBlock, defHS, HS;
  int father, partner, junction, vald0, valr0, childblocks, minblock, recTract, end, beg;
  bool success;
  chrom * chr;
  float p;
  double num;

  // PARENT-PICKING
  father = ancestry[i][t];
  // RECOMBINATION
  if (recombimatrix[i][t] == true) {
    // WHC: childblocks will depends on the position of crossover
    //    childblocks = pointer[prev][father]->b;
    if (ancestry[i][t] % 2 == 0) {
      partner = ancestry[i][t] + 1;
    } else {
      partner = ancestry[i][t] - 1;
    }
    
    chr = pointer[pres][i];

    minblock = minim(pointer[prev][father]->b, pointer[prev][partner]->b);

    if (minblock == 6) {	// WHC: both have 5 blocks, normal crossover
      // here, either from partner or father does not matter
      childblocks = pointer[prev][partner]->b;
      chr->b = childblocks;
    

      //DETERMINE THE HOTSPOT WHERE CROSSOVER WILL OCCUR

      // WHC: Should always consider crossover after deletion happened
      // WHC: probabily will just write a parentpicking() specifically for that phase...
      // WHC: otherwise I think this is workable for now; when a block that doesn't exist happen, just no crossover and copy parental chrom
      p = rand() / ((float) RAND_MAX + 1);
      defHS = numCrossRegions-1;
      num = 1;
      HS = numCrossRegions-1;
      success = 0;
      while(HS >= 0 && success==0){
	num = num-crossRatio[HS];
	if (p < num){
	  defHS = HS-1;
	}
	else {
	  success=1;
	}
	HS--;
      }
      end = crossEnd[defHS];
      beg = crossBegin[defHS];
      if (((minblock * BLOCKLENGTH) < end) && ((minblock * BLOCKLENGTH) > beg)) {
	end = minblock * BLOCKLENGTH;
      }
      if (((minblock * BLOCKLENGTH) >= end) && ((minblock * BLOCKLENGTH) > beg)) {
	// If neither father and partner have mutations in any block
	/*      if ((pointer[prev][father]->mpb[0] == 0) && (pointer[prev][father]->mpb[1] == 0) &&
		(pointer[prev][father]->mpb[2] == 0) && (pointer[prev][partner]->mpb[0] == 0) &&
		(pointer[prev][partner]->mpb[1] == 0) && (pointer[prev][partner]->mpb[2] == 0)) {
	*/

	if ((pointer[prev][father]->mpb[0] == 0) && (pointer[prev][father]->mpb[1] == 0) &&
	    (pointer[prev][father]->mpb[2] == 0) && (pointer[prev][father]->mpb[3] == 0) &&
	    (pointer[prev][father]->mpb[4] == 0) && (pointer[prev][father]->mpb[5] == 0) &&
	    (pointer[prev][partner]->mpb[0] == 0) && (pointer[prev][partner]->mpb[1] == 0) &&
	    (pointer[prev][partner]->mpb[2] == 0) && (pointer[prev][partner]->mpb[3] == 0) &&
	    (pointer[prev][partner]->mpb[4] == 0) && (pointer[prev][partner]->mpb[5] == 0)) {

	  for (j = 0; j < chr->b; j++) {
	    chr->mpb[j] = 0;
	  }
	}// If at least one of them have one or more mutations in one or more blocks
	else {

	  recTract = end - beg;
	  junction = (int) (rand() % (recTract));
	  junction += beg;
	  junctionBlock = (int) (junction / BLOCKLENGTH); // block where junction fell
	  //	   cout << junctionBlock << " ";
	  // WHC: is the number of mutations before junction point
	  valr0 = location(junction - junctionBlock*BLOCKLENGTH, prev, partner, junctionBlock); // junction location in the block "junctionBlock" of partner
	  vald0 = location(junction - junctionBlock*BLOCKLENGTH, prev, father, junctionBlock); // junction location in the block "junctionBlock" of father

	  //COPY THE INFO FROM PARTNER (BLOCK 0)
	  for (j = 0; j < junctionBlock; j++) {
	    chr->mpb[j] = pointer[prev][partner]->mpb[j];
	    for (k = 0; k < pointer[prev][partner]->mpb[j]; k++) {
	      chr->mutation[j][k] = pointer[prev][partner]->mutation[j][k];
	    }
	  }

	  //COPY INFO IN JUNCTION BLOCK (BLOCK 1)
	  //FROM PARTNER

	  // WHC: copy mutation from [partner] until junction point
	  for (k = 0; k < valr0; k++) {
	    chr->mutation[junctionBlock][k] = pointer[prev][partner]->mutation[junctionBlock][k];
	  }
	  //FROM FATHER
	
	  // WHC: copy mutation from [father] from junction point on
	  for (k = 0; k < (pointer[prev][father]->mpb[junctionBlock] - vald0); k++) {
	    chr->mutation[junctionBlock][k + valr0] = pointer[prev][father]->mutation[junctionBlock][k + vald0];
	  }
	  //ESTABLISH MPB
	  chr->mpb[junctionBlock] = valr0 + (pointer[prev][father]->mpb[junctionBlock] - vald0);

	  //COPY INFO FROM BLOCK 2 IF PRESENT

	  // WHC: WRONG??? forgot copying some information???

	
	  // WHC: as from junction point on, all mutation info comes from father...
	  // WHC: (also, previous step only copy mutation info with the junctionBlock!!!
	  for (j = junctionBlock + 1; j < chr->b; j++) {
	    // WHC: as minblock == 6, this should make chr->b = 6, therefore making a for (; j < 6;) loop, copying everything
	    if (j < pointer[prev][father]->b) { //j counts 0,1,2 but .b counts 1,2,3
	      chr->mpb[j] = pointer[prev][father]->mpb[j];
	      for (k = 0 ; k < pointer[prev][father]->mpb[j] ; k++) {
		chr->mutation[j][k] = pointer[prev][father]->mutation[j][k];
	      }
	    } else {
	      chr->mpb[j] = 0;
	    }
	  }
	}
      } else if (((minblock * BLOCKLENGTH) < end) && ((minblock * BLOCKLENGTH) <= beg)) {
	// WHC: this steps looks like a "lazy" step; (NO! I was wrong! Imagine if a hot block on block 3 with 100% proposion!
	// WHC: (then, if we don't have block 3, we don't have recombination!)

	// WHC: WRONG!! copychr is wrong!!!
	/* ================================================================ */
	//	copychr(prev, father, pres, i);
	copychr_for_phaseVI(prev, father, pres, i);
      }
    } else if (minblock == 5) {
      // WHC: this is when one chrom has 4 blocks while another has 5 blocks; tricky one

      // WHC: let's decide HS0(block 0, 1, 2), HS1(block 3) and HS2(block 4), which one will have
      p = rand() / ((float) RAND_MAX + 1);
      defHS = numCrossRegions-1; // as numCrossRegions == numHS == 3, this defHS = 2 here
      num = 1;
      HS = numCrossRegions-1;
      success = 0;
      while(HS >= 0 && success==0){
	// WHC: so, defHS is the final index(0, 1, 2) of crossoverRatio etc.
	num = num-crossRatio[HS];
	if (p < num){
	  defHS = HS-1;
	}
	else {
	  success=1;
	}
	HS--;
      }

      //      if (defHS == 0) {		// WHC: means crossover will happen on block 0 - 2, which needs to escape block 2, as it is absent
	// WHC: still, partner decides # of child's block
	// WHC: depends on father here! on partner in else if clause

	/*************************************************************************************************************************
        WHC: as stated before, when defHS == 0, information before junction point comes from partner
        chr->b = #blocks for father
	**************************************************************************************************************************/

      if (defHS == 0 || defHS == 1) {		// WHC: means crossover happens on block 0 or 1; or block 2
	childblocks = pointer[prev][father]->b;
	chr->b = childblocks;

	// WHC: crossover could only happen on block 0 or 1
	//	end = 3 * BLOCKLENGTH;
	end = crossEnd[defHS];
	beg = crossBegin[defHS]; // WHC: this is to say, beg = 0
	if (beg > 2 * BLOCKLENGTH) { cout << "this could not be!"; exit(0); }
	if (end > 3 * BLOCKLENGTH) { cout << "this could not be like this!"; exit(0); }

	// WHC: crossover may happen in either block 0 or block 1; but do NOT know how many blocks father has; so has to test all

	if ((pointer[prev][father]->mpb[0] == 0) && (pointer[prev][father]->mpb[1] == 0) &&
	    (pointer[prev][father]->mpb[2] == 0) && (pointer[prev][father]->mpb[3] == 0) &&
	    (pointer[prev][father]->mpb[4] == 0) && (pointer[prev][father]->mpb[5] == 0) &&
	    (pointer[prev][partner]->mpb[0] == 0) && (pointer[prev][partner]->mpb[1] == 0) &&
	    (pointer[prev][partner]->mpb[2] == 0) && (pointer[prev][partner]->mpb[3] == 0) &&
	    (pointer[prev][partner]->mpb[4] == 0) && (pointer[prev][partner]->mpb[5] == 0)) {
	  for (j = 0; j < chr->b; j++) {
	    chr->mpb[j] = 0;
	  }
	}// If at least one of them have one or more mutations in one or more blocks
	else {

	  recTract = end - beg;
	  junction = (int) (rand() % (recTract));
	  junction += beg;
	  junctionBlock = (int) (junction / BLOCKLENGTH); // block where junction fell
	  if (junctionBlock > 2) { cout << "this shouldn't happne too!\n"; exit(0); } // WHC: just a doulbe-check
	  
	  //	   cout << junctionBlock << " ";
	  // WHC: is the number of mutations before junction point
	  valr0 = location(junction - junctionBlock*BLOCKLENGTH, prev, partner, junctionBlock); // junction location in the block "junctionBlock" of partner
	  vald0 = location(junction - junctionBlock*BLOCKLENGTH, prev, father, junctionBlock); // junction location in the block "junctionBlock" of father

	  //COPY THE INFO FROM PARTNER (BLOCK 0)
	  for (j = 0; j < junctionBlock; j++) {
	    // WHC: if partner has dup_1, then it is fine
	    // WHC: else if it doesn't have dup_1, them mpb[2] should be 0, so this should be fine; need double-check
	    chr->mpb[j] = pointer[prev][partner]->mpb[j];
	    for (k = 0; k < pointer[prev][partner]->mpb[j]; k++) {
	      chr->mutation[j][k] = pointer[prev][partner]->mutation[j][k];
	    }
	  }

	  // WHC: copy mutation from [partner] until junction point
	  for (k = 0; k < valr0; k++) {
	    chr->mutation[junctionBlock][k] = pointer[prev][partner]->mutation[junctionBlock][k];
	  }
	  //FROM FATHER
	
	  // WHC: copy mutation from [father] from junction point on; until the end of junctionBlock
	  for (k = 0; k < (pointer[prev][father]->mpb[junctionBlock] - vald0); k++) {
	    chr->mutation[junctionBlock][k + valr0] = pointer[prev][father]->mutation[junctionBlock][k + vald0];
	  }
	  //ESTABLISH MPB
	  chr->mpb[junctionBlock] = valr0 + (pointer[prev][father]->mpb[junctionBlock] - vald0);


	  // WHC: WRONG??? forgot copying some information???

	
	  // WHC: as from junction point on, all mutation info comes from father...
	  // WHC: (also, previous step only copy mutation info with the junctionBlock!!!
	  //	  for (j = junctionBlock + 1; j < chr->b; j++) {
	  // WHC: here, should change; as chr->b may be 4 while #blocks for fathter = 5
	  // j should == junctionBlock now
	  j = junctionBlock + 1;

	  if ((j == 1) && (pointer[prev][father]->b == 6)) { //j counts 0,1,2 but .b counts 1,2,3
	    // WHC: WRONG at first!! should copy father's information in ALL block 1, 2, 3, and 4!!!
	    for (int temp = 1; temp < 6; ++temp) {
	      chr->mpb[temp] = pointer[prev][father]->mpb[temp];
	      for (k = 0 ; k < pointer[prev][father]->mpb[temp] ; k++) {
		      chr->mutation[temp][k] = pointer[prev][father]->mutation[temp][k];
	      }
	    }
	    /* WHC: wrong at the first time, only copied information in block 1 & 2, not in block 3, 4
	    chr->mpb[1] = pointer[prev][father]->mpb[1];
	    chr->mpb[2] = pointer[prev][father]->mpb[2];
	    for (k = 0 ; k < pointer[prev][father]->mpb[1] ; k++) {
	      chr->mutation[1][k] = pointer[prev][father]->mutation[1][k];
	    }
	    for (k = 0 ; k < pointer[prev][father]->mpb[2] ; k++) {
	      chr->mutation[2][k] = pointer[prev][father]->mutation[2][k];
	    }
	    */
	  } else if (j == 1 && pointer[prev][father]->b == 5) {

	    for (int temp = 1; temp < 6; ++temp) {
	      chr->mpb[temp] = pointer[prev][father]->mpb[temp];
	      for (k = 0 ; k < pointer[prev][father]->mpb[temp] ; k++) {
		//		if (temp == 2) { continue; } , don't need this, as pointer[][]->mpb[2] == 0 already
		chr->mutation[temp][k] = pointer[prev][father]->mutation[temp][k];
	      }
	    }

	    /*
	    chr->mpb[1] = pointer[prev][father]->mpb[1];
	    // WHC: important!!! need also assign pointer[prev][father]->mpb[2] (which is 0) to child!!!!!!!! for viod statistics_for_phaseVI()
	    chr->mpb[2] = pointer[prev][father]->mpb[2];
	    for (k = 0 ; k < pointer[prev][father]->mpb[1] ; k++) {
	      chr->mutation[1][k] = pointer[prev][father]->mutation[1][k];
	    }
	    */
	    
	  } else if (j == 2) {
	    // WHC: do nothing
	    // WCH: NO!!! need to copy information in block 2, 3, 4 to child chrom!! father->mpb[2] = 0 need to be copied too!
	    for (int temp = 2; temp < 6; ++temp) {
	      chr->mpb[temp] = pointer[prev][father]->mpb[temp];
	      for (k = 0 ; k < pointer[prev][father]->mpb[temp] ; k++) {
		// WHC: don't skip block 2 this time
		chr->mutation[temp][k] = pointer[prev][father]->mutation[temp][k];
	      }
	    }
	    
	  } else if (j == 3) {
	    for (int temp = 3; temp < 6; ++temp) {
	      chr->mpb[temp] = pointer[prev][father]->mpb[temp];
	      for (k = 0 ; k < pointer[prev][father]->mpb[temp] ; k++) {
		// WHC: don't skip block 2 this time
		chr->mutation[temp][k] = pointer[prev][father]->mutation[temp][k];
	      }
	    }
	  } else {
	    cout << "impossible!\n";
	    exit(0);
	    //chr->mpb[j] = 0;
	  }
	}
	//      } else if (defHS == 1 || defHS == 2) { // means crossover will happen on block 3 - 4, needs to pick #blocks for child chrom
      } else if (defHS == 3 || defHS == 4) {
	//      if (defHS == 1 || defHS == 2) { // means crossover will happen on block 3 - 4, needs to pick #blocks for child chrom
	// WHC: I incoporated defHS == 0 here, just need to adjust beg and end
	// WHC: cannot merge them; because if crossover happen before block 2, chr->b depends on father; otherwise on partner
	
	// WHC: as crossover happens after dup_1, #blocks(child) = #block(partner); as designed in parentpicking(), left side comes from partner, while right side comes from father


	/*
	  if crossover happens in block 3 or block 4, all information before junction point comes from father,
	  to make sure chr->b = pointer[prev][father]->b
	*/
	
	//	childblocks = pointer[prev][partner]->b;
	childblocks = pointer[prev][father]->b;
	chr->b = childblocks;
    

      //DETERMINE THE HOTSPOT WHERE CROSSOVER WILL OCCUR

      // WHC: Should always consider crossover after deletion happened
      // WHC: probabily will just write a parentpicking() specifically for that phase...
      // WHC: otherwise I think this is workable for now; when a block that doesn't exist happen, just no crossover and copy parental chrom
      // as defHS already decided
      if (defHS == 0 || defHS == 1) {
	cout << "this should not happen.\n";
	exit(0);
	end = 3 * BLOCKLENGTH;
	beg = crossBegin[defHS];
	if (beg != 0) { cout << "this could not be right.\n"; exit(0); }
      } else {
	end = crossEnd[defHS];
	beg = crossBegin[defHS];
      }
      // WHC: the following if clause no longer fits here
      //      if (((minblock * BLOCKLENGTH) < end) && ((minblock * BLOCKLENGTH) > beg)) {
      //	end = minblock * BLOCKLENGTH;
      //      }
      

	// If neither father and partner have mutations in any block
	/*      if ((pointer[prev][father]->mpb[0] == 0) && (pointer[prev][father]->mpb[1] == 0) &&
		(pointer[prev][father]->mpb[2] == 0) && (pointer[prev][partner]->mpb[0] == 0) &&
		(pointer[prev][partner]->mpb[1] == 0) && (pointer[prev][partner]->mpb[2] == 0)) {
	*/
      // WHC: as when lose_of_duplication(), I set mpb[0] == 0, this should be fine???
	if ((pointer[prev][father]->mpb[0] == 0) && (pointer[prev][father]->mpb[1] == 0) &&
	    (pointer[prev][father]->mpb[2] == 0) && (pointer[prev][father]->mpb[3] == 0) &&
	    (pointer[prev][father]->mpb[4] == 0) && (pointer[prev][father]->mpb[5] == 0) &&
	    (pointer[prev][partner]->mpb[0] == 0) && (pointer[prev][partner]->mpb[1] == 0) &&
	    (pointer[prev][partner]->mpb[2] == 0) && (pointer[prev][partner]->mpb[3] == 0) &&
	    (pointer[prev][partner]->mpb[4] == 0) && (pointer[prev][partner]->mpb[5] == 0)) {
	  for (j = 0; j < chr->b; j++) {
	    chr->mpb[j] = 0;
	  }
	}// If at least one of them have one or more mutations in one or more blocks
	else {

	  recTract = end - beg;
	  junction = (int) (rand() % (recTract));
	  junction += beg;
	  junctionBlock = (int) (junction / BLOCKLENGTH); // block where junction fell
	  if (defHS > 0 && junctionBlock < 4) { cout << "this shouldn't happne!\n"; exit(0); } // WHC: just a doulbe-check
	  
	  //	   cout << junctionBlock << " ";
	  // WHC: is the number of mutations before junction point
	  valr0 = location(junction - junctionBlock*BLOCKLENGTH, prev, partner, junctionBlock); // junction location in the block "junctionBlock" of partner
	  vald0 = location(junction - junctionBlock*BLOCKLENGTH, prev, father, junctionBlock); // junction location in the block "junctionBlock" of father

	  //COPY THE INFO FROM PARTNER (BLOCK 0)
	  // WHC: NO, from father
	  for (j = 0; j < junctionBlock; j++) {
	    // WHC: if partner has dup_1, then it is fine
	    // WHC: else if it doesn't have dup_1, them mpb[2] should be 0, so this should be fine; need double-check
	    chr->mpb[j] = pointer[prev][father]->mpb[j];
	    for (k = 0; k < pointer[prev][father]->mpb[j]; k++) {
	      chr->mutation[j][k] = pointer[prev][father]->mutation[j][k];
	    }
	  }


	  // WHC: copy mutation from [FATHER] until junction point
	    //	  for (k = 0; k < valr0; k++) {
	    for (k = 0; k < vald0; k++) {
	    chr->mutation[junctionBlock][k] = pointer[prev][father]->mutation[junctionBlock][k];
	    }
	  //FROM FATHER
	
	  // WHC: copy mutation from [PARTNER] from junction point on
	    for (k = 0; k < (pointer[prev][partner]->mpb[junctionBlock] - valr0); k++) {
	      chr->mutation[junctionBlock][k + vald0] = pointer[prev][partner]->mutation[junctionBlock][k + valr0];
	    }
	  //ESTABLISH MPB
	  chr->mpb[junctionBlock] = vald0 + (pointer[prev][partner]->mpb[junctionBlock] - valr0);

	  //COPY INFO FROM BLOCK 2 IF PRESENT

	  // WHC: WRONG??? forgot copying some information???

	
	  // WHC: as from junction point on, all mutation info comes from father...
	  // WHC: (also, previous step only copy mutation info with the junctionBlock!!!
	  //	  for (j = junctionBlock + 1; j < chr->b; j++) {
	  // WHC: here, should change; as chr->b may be 4 while #blocks for fathter = 5
	  // j should == junctionBlock now
	  j = junctionBlock + 1;

	  if (j == 5) { //j counts 0,1,2 but .b counts 1,2,3
	    chr->mpb[5] = pointer[prev][partner]->mpb[5];
	    for (k = 0 ; k < pointer[prev][partner]->mpb[5] ; k++) {
	      chr->mutation[5][k] = pointer[prev][partner]->mutation[5][k];
	    }
	  } else if (j == 6) {
	    // WHC: do nothing
	  } else {
	    cout << "impossible!\n";
	    exit(0);
	    //chr->mpb[j] = 0;
	  }


	}
	
      } else if (defHS == 2) {
	// WHC: this means no crossover
	copychr_for_phaseVI(prev, father, pres, i);
      }
      
    } else {
      cout << "something wrong in parentpicking_for_phaseVI() here.\n";
      exit(0);
    }
  }// NO RECOMBINATION
  else {
    // WHC: WRONG!!! Cannot use original copychr, as it uses pointer[prev][father]->b for a for loop
    //    copychr(prev, father, pres, i);
    copychr_for_phaseVI(prev, father, pres, i);
  }

}

// WHC: the end of parentpicking_for_phaseVI()
/* ================================================================ */
void duplication(int i,int prev, bool from) {
  int k;
  int tempMutCount = 0;
  int adam = i;

  if (from == 0) {
    if (i % 2 == 0) {
      adam = i + 1;
    } else {
      adam = i - 1;
    }
  }
  // WHC: pointer[][] contains pointers to table[]; and table[] contains struct chrom;
  // WHC: so pointer[][][i] represents ith chromes, because array name == pointer
  if (pointer[prev][i][0].b == 2) {
    //    pointer[prev][i][0].b++;
    pointer[prev][i][0].b = 4;
  }
  //  pointer[prev][i][0].mpb[2] = pointer[prev][adam][0].mpb[0];
  pointer[prev][i][0].mpb[3] = pointer[prev][adam][0].mpb[0];
  // WHC: spacer info
  pointer[prev][i]->mpb[2] = 0;
  if (pointer[prev][adam]->mpb[2] != 0) {
    cout << "The spacer should have 0 mutations before duplication.\n";
    exit(0);
  }
  
  for (k = 0; k < pointer[prev][adam][0].mpb[0]; k++) {
    pointer[prev][i][0].mutation[3][k] = pointer[prev][adam][0].mutation[0][k];
  }
  for (k = 0; k < MutCount; k++) {
    if (muttable[k].block == 0) {
      muttable[MutCount + tempMutCount].position = muttable[k].position;
      muttable[MutCount + tempMutCount].block = 3;
      tempMutCount++;
    }
  }
  MutCount += tempMutCount;
}

/* ================================================================================ */
// WHC: new duplication_2() for phaseIV()

void duplication_2(int i,int prev, bool from) {
  int k;
  int tempMutCount = 0;
  int adam = i;

  if (from == 0) {
    if (i % 2 == 0) {
      adam = i + 1;
    } else {
      adam = i - 1;
    }
  }
  // WHC: pointer[][] contains pointers to table[]; and table[] contains struct chrom;
  // WHC: so pointer[][][i] represents ith chromes, because array name == pointer
  if (pointer[prev][i]->b == 4) {
    pointer[prev][i][0].b = 6;
  }

  // WHC: also, will assume dup_2 comes from ori now... may need to re-consider
  
  pointer[prev][i][0].mpb[5] = pointer[prev][adam][0].mpb[0];
  pointer[prev][i]->mpb[4] = 0;
  for (k = 0; k < pointer[prev][adam][0].mpb[0]; k++) {
    pointer[prev][i][0].mutation[5][k] = pointer[prev][adam][0].mutation[0][k];
  }

  // WHC: need to understand this before changing; changed, need double-check
  for (k = 0; k < MutCount; k++) {
    // WHC: this may be the causes of thousands of tries!!!! Should also copy dup_1's informations!!! (WRONG!!!)
    
    //    if (muttable[k].block == 0 || muttable[k].block == 2) {
    if (muttable[k].block == 0) {
      /*
	WHC: WAIT!! But does I copy twice?? If block 0 and block 2 share all mutation information, am I copying twice???
	WHC: IMPORTANT!!!!
	WHC: after changing this, the redundency in muttable[].block = 4 is eliminated!!!
       */
      
      muttable[MutCount + tempMutCount].position = muttable[k].position;
      muttable[MutCount + tempMutCount].block = 5;
      tempMutCount++;
    }
  }
  MutCount += tempMutCount;
}

// WHC: new duplication_2() for phaseIV()
/* ================================================================================ */

/* ================================================================================ */
// WHC: new lose_of_duplication for phaseVI()

void lose_of_duplication(int eva, int prev) {
  // WHC: need to eliminate muttable[][]'s information that came from this block????
  // WHC: as even void mutation() only added new mutation into muttable, and lost mutations not erased until FSL()
  // WHC: I think it is fine to just ignore this
  // WHC: also, should consider making mpb[2] == 0?
  
  if (pointer[prev][eva]->b == 6) {
    // WHC: actually chr->b is just a indicator now. still has 6 blocks, and mpb[3] = 0 (dup_1); that's why should use mutation_for_phaseVI()
    pointer[prev][eva]->b = 5;

    // WHC: this step is VERY important, as many problems have been caused by this!!!
    pointer[prev][eva]->mpb[3] = 0;
    
  } else {
    cout << "something wrong in lose_of_duplication().\n";
    exit(0);
  }
}

// WHC: new lose_of_duplication for phaseVI()
/* ================================================================================ */

void mutation(float probability, int i, int pres) {
  int j, k, position, val, mutEvents;
  float p;
  chrom * chr;
  chr = pointer[pres][i];


  /* WHC: For test purpose; under realistic scenario, multihit[] == true should not be larger than BLOCKLENGTH
  int count_multihit = 0, count_multihit_single = 0;
  for (int t = 0; t < BLOCKLENGTH; ++t) {
    if (multihit[t] == true) {
      ++count_multihit;
    }
    if (multihit_single[t] == true) {
      ++count_multihit_single;
    }
  }
  if (count_multihit >= BLOCKLENGTH || count_multihit_single >= BLOCKLENGTH) {
    cout << "The multihit or multihit_single is larger than BLOCKLENGTH.\nWill have an infinate loop!\n";
    exit(0);
  }
  */
  
  for (j = 0; j < chr->b; j++) {

    if (j == 2 || j == 4) {
      continue; 	// WHC: skip block 2 or 4, the two spacers
    }
    
    mutEvents=0;
    p = rand() / ((float) RAND_MAX + 1);
    if (p < (probability * BLOCKLENGTH)) { mutEvents++;
      if (p < pow((probability * BLOCKLENGTH), 2)) { mutEvents++;
	if (p < pow((probability * BLOCKLENGTH), 3)) { mutEvents++;
	  if (p < pow((probability * BLOCKLENGTH),4)) { mutEvents++;
	  }
	}
      }
    }
    if (j != 1) {

      /*
      if (duFreq_2 == true) {
	how_many_new_mutations_in_dup += mutEvents;
      }
      */
      
      for(int muts=0; muts < mutEvents; muts++){
	// Randomly choose one NEW mutation position
	do {
	  position = (int) (rand() % (BLOCKLENGTH));
	} while (multihit[position] == true);
	multihit[position] = true;


	// If there are no previous mutations or if the new mutation fall behind all previous ones...,
	if (chr->mpb[j] == 0 || position > chr->mutation[j][chr->mpb[j] - 1]) {
	  // Simply add the mutation
	  chr->mutation[j][chr->mpb[j]] = position;
	  chr->mpb[j]++;
	}// Otherwise, shift the mutations one position so that mutations appear ordered in the chr.mutation array
	// WHC: otherwise shift many positions; I think could be done by stable_sort()?
	else {
	  val = location(position, pres, i, j);
	  for (k = chr->mpb[j]; k > val; k--) {
	    // WHC: from last mutation, move everything 1 position backword, until new mutation
	    chr->mutation[j][k] = chr->mutation[j][k - 1];
	  }
	  chr->mutation[j][val] = position;
	  chr->mpb[j]++;
	}
	// ELABORATES MUTTABLE
	muttable[MutCount].position = position;
	muttable[MutCount].block = j;
	MutCount++;


	// WHC: I think this step is for making sure new mutation could only arise on non-polymorphic positions, on both blocks
	// IF MUTATION APPEARS IN BLOCK 0 OR 2, COPY THE MUTATION TO THE OTHER BLOCK
	// if (j == 0 && duFreq == true) {
	if (j == 0 && duFreq == true && duFreq_2 == false) {
	  // WHC: in phaseII() but not in phaseIV()
	  muttable[MutCount].block = 3;
	  muttable[MutCount].position = muttable[MutCount - 1].position;
	  MutCount++;
	} else if (j == 0 && duFreq_2 == true) {
	  // WHC: if is in phaseIV(), then should copy information from 0 to 2 and 4
	  // WHC: copy to dup_1
	  muttable[MutCount].block = 3;
	  muttable[MutCount].position = muttable[MutCount - 1].position;
	  MutCount++;

	  // WHC: copy to dup_2
	  muttable[MutCount].block = 5;
	  muttable[MutCount].position = muttable[MutCount - 2].position;
	  MutCount++;
	}
      
	//      if (j == 2) {
	if (j == 3 && duFreq_2 == false) {
	  // WHC: in phaseII() but not in phaseIV()
	  muttable[MutCount].block = 0;
	  muttable[MutCount].position = muttable[MutCount - 1].position;
	  MutCount++;
	} else if (j == 3 && duFreq_2 == true) {
	  // WHC: in phaseIV()
	  // WHC: or in phaseV()
	  // WHC: copy to dup_1
	  muttable[MutCount].block = 0;
	  muttable[MutCount].position = muttable[MutCount - 1].position;
	  MutCount++;

	  // WHC: copy to dup_2
	  muttable[MutCount].block = 5;
	  muttable[MutCount].position = muttable[MutCount - 2].position;
	  MutCount++;
	}

	// WHC: FORGOT the j == 4!!!! Causing trouble in muttable[]!!!! BUT, I am so happy figuring this out!!!

	if (j == 5) {
	  // WHC: must be in or after phaseIV()
	  muttable[MutCount].block = 0;
	  muttable[MutCount].position = muttable[MutCount - 1].position;
	  MutCount++;

	  // WHC: copy to dup_2
	  muttable[MutCount].block = 3;
	  muttable[MutCount].position = muttable[MutCount - 2].position;
	  MutCount++;
	}
      
      }

	

    } else if (j == 1) {
      for(int muts=0; muts < mutEvents; muts++){
	// Randomly choose one NEW mutation position
	do {
	  position = (int) (rand() % (BLOCKLENGTH));
	} while (multihit_single[position] == true);
	multihit_single[position] = true;


	// If there are no previous mutations or if the new mutation fall behind all previous ones...,
	if (chr->mpb[j] == 0 || position > chr->mutation[j][chr->mpb[j] - 1]) {
	  // Simply add the mutation
	  chr->mutation[j][chr->mpb[j]] = position;
	  chr->mpb[j]++;
	}// Otherwise, shift the mutations one position so that mutations appear ordered in the chr.mutation array
	// WHC: otherwise shift many positions; I think could be done by stable_sort()?
	else {
	  val = location(position, pres, i, j);
	  for (k = chr->mpb[j]; k > val; k--) {
	    // WHC: from last mutation, move everything 1 position backword, until new mutation
	    chr->mutation[j][k] = chr->mutation[j][k - 1];
	  }
	  chr->mutation[j][val] = position;
	  chr->mpb[j]++;
	}
	// ELABORATES MUTTABLE
	muttable[MutCount].position = position;
	muttable[MutCount].block = j;
	MutCount++;

      }

    } else {
      cout << "WRONG\n";
      exit(0);
    }

  }

}


/* ================================================================ */
// WHC: new mutation() for phaseVI()

void mutation_for_phaseVI(float probability, int i, int pres) {
  int j, k, position, val, mutEvents;
  float p;
  chrom * chr;
  chr = pointer[pres][i];

  // WHC: just slightly adjust this void mutation(), by skiping j == 2 if chr->b == 4
  if (loseFreq == false) { cout << "This mutation_for_phaseVI() is only used for phaseVI().\n"; exit(0); }
  
  for (j = 0; j < 6; j++) {
    if (chr->b == 5 && j == 3) { continue; }
    // WHC: just skip block 2 and 4, the two spacers
    if (j == 2 || j == 4) { continue; }

    mutEvents=0;
    p = rand() / ((float) RAND_MAX + 1);
    if (p < (probability * BLOCKLENGTH)) { mutEvents++;
      if (p < pow((probability * BLOCKLENGTH), 2)) { mutEvents++;
	if (p < pow((probability * BLOCKLENGTH), 3)) { mutEvents++;
	  if (p < pow((probability * BLOCKLENGTH),4)) { mutEvents++;
	  }
	}
      }
    }

    if (j != 1) {
      
      //      how_many_new_mutations_in_dup += mutEvents;
      
    for(int muts=0; muts < mutEvents; muts++){
      // Randomly choose one NEW mutation position
      do {
	position = (int) (rand() % (BLOCKLENGTH));
      } while (multihit[position] == true);
      multihit[position] = true;

      // If there are no previous mutations or if the new mutation fall behind all previous ones...,
      if (chr->mpb[j] == 0 || position > chr->mutation[j][chr->mpb[j] - 1]) {
	// Simply add the mutation
	chr->mutation[j][chr->mpb[j]] = position;
	chr->mpb[j]++;
      }// Otherwise, shift the mutations one position so that mutations appear ordered in the chr.mutation array
      // WHC: otherwise shift many positions; I think could be done by stable_sort()?
      else {
	val = location(position, pres, i, j);
	for (k = chr->mpb[j]; k > val; k--) {
	  // WHC: from last mutation, move everything 1 position backword, until new mutation
	  chr->mutation[j][k] = chr->mutation[j][k - 1];
	}
	chr->mutation[j][val] = position;
	chr->mpb[j]++;
      }
      // ELABORATES MUTTABLE
      muttable[MutCount].position = position;
      muttable[MutCount].block = j;
      MutCount++;


      // WHC: I think this step is for making sure new mutation could only arise on non-polymorphic positions, on both blocks
      // IF MUTATION APPEARS IN BLOCK 0 OR 2, COPY THE MUTATION TO THE OTHER BLOCK
      // if (j == 0 && duFreq == true) {
      if (j == 0 && duFreq == true && duFreq_2 == false) {
	// WHC: in phaseII() but not in phaseIV()
	muttable[MutCount].block = 3;
	muttable[MutCount].position = muttable[MutCount - 1].position;
	MutCount++;
      } else if (j == 0 && duFreq_2 == true) {
	// WHC: if is in phaseIV(), then should copy information from 0 to 2 and 4
	// WHC: copy to dup_1
	muttable[MutCount].block = 3;
	muttable[MutCount].position = muttable[MutCount - 1].position;
	MutCount++;

	// WHC: copy to dup_2
	muttable[MutCount].block = 5;
	muttable[MutCount].position = muttable[MutCount - 2].position;
	MutCount++;
      }
      
      //      if (j == 2) {
      if (j == 3 && duFreq_2 == false) {
	// WHC: in phaseII() but not in phaseIV()
	muttable[MutCount].block = 0;
	muttable[MutCount].position = muttable[MutCount - 1].position;
	MutCount++;
      } else if (j == 3 && duFreq_2 == true) {
	// WHC: in phaseIV()
	// WHC: or in phaseV()
	// WHC: copy to dup_1
	muttable[MutCount].block = 0;
	muttable[MutCount].position = muttable[MutCount - 1].position;
	MutCount++;

	// WHC: copy to dup_2
	muttable[MutCount].block = 5;
	muttable[MutCount].position = muttable[MutCount - 2].position;
	MutCount++;
      }

      // WHC: FORGOT the j == 4!!!! Causing trouble in muttable[]!!!! BUT, I am so happy figuring this out!!!

      if (j == 5) {
	// WHC: must be in or after phaseIV()
	muttable[MutCount].block = 0;
	muttable[MutCount].position = muttable[MutCount - 1].position;
	MutCount++;

	// WHC: copy to dup_2
	muttable[MutCount].block = 3;
	muttable[MutCount].position = muttable[MutCount - 2].position;
	MutCount++;
      }
      
    }
    } else if (j == 1) {
      for(int muts=0; muts < mutEvents; muts++){
	// Randomly choose one NEW mutation position
	do {
	  position = (int) (rand() % (BLOCKLENGTH));
	} while (multihit_single[position] == true);
	multihit_single[position] = true;

	// If there are no previous mutations or if the new mutation fall behind all previous ones...,
	if (chr->mpb[j] == 0 || position > chr->mutation[j][chr->mpb[j] - 1]) {
	  // Simply add the mutation
	  chr->mutation[j][chr->mpb[j]] = position;
	  chr->mpb[j]++;
	}// Otherwise, shift the mutations one position so that mutations appear ordered in the chr.mutation array
	// WHC: otherwise shift many positions; I think could be done by stable_sort()?
	else {
	  val = location(position, pres, i, j);
	  for (k = chr->mpb[j]; k > val; k--) {
	    // WHC: from last mutation, move everything 1 position backword, until new mutation
	    chr->mutation[j][k] = chr->mutation[j][k - 1];
	  }
	  chr->mutation[j][val] = position;
	  chr->mpb[j]++;
	}
	// ELABORATES MUTTABLE
	muttable[MutCount].position = position;
	muttable[MutCount].block = j;
	MutCount++;
      }
    } else {
      cout << "wrong\n";
      exit(0);
    }
  }

}

/* ================================================================ */

// void conversion(float probability, int t, int i, int pres, float donorRatio, float sameDifIGC) {

// WHC: add a new step here, to select a pair of duplicated genes, i.e. ori & dup_1 or dup_1 & dup_2 etc. to do conversion
void conversion(float probability, int t, int i, int pres, float (*p_donorRatio)[5], float sameDifIGC) {
  int k, partner, donor, receptor, chrDonor, chrReceptor, IGC, vald0, valdf, valr0, valrf, val, junction, tractlength, eqmut, otherk, differences;

  float p;
  chrom * chr1, * chr2;
  chr1 = pointer[pres][i];
  IGC = 0;



  
  // If the chr has duplicated block
  //  if (chr1[0].b == 3) {
  // WHC: as long as #blocks >= 3
  // WHC: REMEMBER, #blocks also determine which 2 pairs of duplications could be selected
  
  if (chr1->b >= 4) {
    if((sameDifIGC!=1) && (IGCmatrix[i][t]==true)){ // IGC on different chroms
      IGC = 1;

      if (i % 2 == 0) {
	partner = (i + 1);
      } else {
	partner = (i - 1);
      }			
      chr2 = pointer[pres][partner];

      /* ================================================================================================ */
      // beginning of new donor/receptor selection method
      
      // WHC: randomly choose a pair of duplicated blocks, depending on the number of blocks chrom has
      // WHC: also, return value of donorRatio
      int block_1, block_2;
      float donorRatio;
      
      if (chr1->b == 4 && chr2->b <= 4) {		// WHC: pick 2 pairs
	block_1 = ori_index;
	block_2 = dup_1;
	donorRatio = p_donorRatio[0][2];
      } else if (chr1->b == 6 || chr2->b == 6) {	// WHC: could pick ori, dup_1, dup_2 blocks; means already in phaseIV(), at least 3 blocks for any chrom
	// WHC: randomly generate 0, 1, 2
	// WHC: which correspond to ori&dup_1, ori&dup_2 and dup_1&dup_2 pairs

	// WHC: pick_a_pair() is a function, that picks 0, 1 or 2, with desired propotions
	int which_pair = pick_a_pair(1);
	//	int which_pair = rand() % 3;

	if (which_pair == 0) {	// WHC: ori&dup_1 pair
	  block_1 = ori_index;
	  block_2 = dup_1;
	  donorRatio = p_donorRatio[0][2];
	} else if (which_pair == 1) { // WHC: ori&dup_2 pair
	  block_1 = ori_index;
	  block_2 = dup_2;
	  donorRatio = p_donorRatio[0][4];
	} else if (which_pair == 2) { // WHC: dup_1&dup_2 pair
	  block_1 = dup_1;
	  block_2 = dup_2;
	  donorRatio = p_donorRatio[2][4];
	} else {			// WHC: wrong
	  cout << "something is wrong in void conversion()\n";
	  exit(0);
	}
      } else {
	cout << "wrong on void conversion() 2\n";
	exit(0);
      }

      // Determines which block will be the donor and which will be the receptor
      p = rand() / ((float) RAND_MAX + 1);
      if (p < donorRatio) {
	donor = block_1;
	receptor = block_2;
      } else {
	donor = block_2;
	receptor = block_1;
      }
      
      if (chr2->b == 6 && chr1->b == 6) {
	p = rand() / ((float) RAND_MAX + 1); // WHC: randomly choose from which chromosome to which
	if (p < 0.5) {
	  chrDonor = i;
	  chrReceptor = partner;
	} else {
	  chrDonor = partner;
	  chrReceptor = i;
	}
      } else if (chr2->b == 6 && chr1->b == 4) {
	if (block_2 == dup_2) {	// only chr2 can be donor
	  chrDonor = partner;
	  chrReceptor = i;
	} else {		// randomly chosen; WHC: as block_2 = dup_1... (WHC: wise to make block_2 > block_1)
	  p = rand() / ((float) RAND_MAX + 1); // WHC: randomly choose from which chromosome to which
	  if (p < 0.5) {
	    chrDonor = i;
	    chrReceptor = partner;
	  } else {
	    chrDonor = partner;
	    chrReceptor = i;
	  }
	}
      } else if (chr2->b == 4) { // WHC: chr2 only contains ori and dup_1 (phase which loses dup_1 will not be considered here
	// chr1->b could either be 3 or 5
	if (block_2 == dup_2) {
	  chrDonor = i;
	  chrReceptor = partner;
	} else {
	  p = rand() / ((float) RAND_MAX + 1); // WHC: randomly choose from which chromosome to which
	  if (p < 0.5) {
	    chrDonor = i;
	    chrReceptor = partner;
	  } else {
	    chrDonor = partner;
	    chrReceptor = i;
	  }
	}
	
      } else if (chr2->b == 2) {
	// WHC: may be a source of errors, careful
	//	if(donor == 2){
	if (donor == dup_1 || donor == dup_2) { // WHC: this chrom does not have duplications
	  chrDonor = i;
	  chrReceptor = partner;
	} else {
	  chrDonor = partner;
	  chrReceptor = i;
	}
      }

      // end of new donor/receptor selection method
      /* ================================================================================================ */
      
    } else {
      // WHC: if IGCmatrix[i][t] == false, this means that no IGC between chroms will happen, due to void genealogy()
      // WHC: but IGC could still happen within chrom
      p = rand() / ((float) RAND_MAX + 1);
      // WHC: should distinguish between 2 duplications and 3 duplications
      if (chr1->b == 4) {
	int block_1 = ori_index;
	int block_2 = dup_1;
	float donorRatio = p_donorRatio[0][2];

	if (p < (2 * probability * BLOCKLENGTH * sameDifIGC)) {
	  IGC = 1;
	  // Determines which block will be the donor and which will be the receptor
	  p = rand() / ((float) RAND_MAX + 1);
	  if (p < donorRatio) {
	    donor = block_1;
	    receptor = block_2;
	  } else {
	    donor = block_2;
	    receptor = block_1;
	  }
	  chrDonor = i;
	  chrReceptor = i;
	}
      } else if (chr1->b == 6) {
	// WHC: I should also add chr1->b == 5 (lost dup_1!!!)
	if (p < (3 * probability * BLOCKLENGTH * sameDifIGC)) { // WHC: should this be 3? Because we have 3 blocks...
	  IGC = 1;
	  // Determines which block will be the donor and which will be the receptor

	  
	  int which_pair = pick_a_pair(1);
	  //	  int which_pair = rand() % 3;
	  int block_1, block_2;
	  float donorRatio;
	  if (which_pair == 0) {	// WHC: ori&dup_1 pair
	    block_1 = ori_index;
	    block_2 = dup_1;
	    donorRatio = p_donorRatio[0][2];
	  } else if (which_pair == 1) { // WHC: ori&dup_2 pair
	    block_1 = ori_index;
	    block_2 = dup_2;
	    donorRatio = p_donorRatio[0][4];
	  } else if (which_pair == 2) { // WHC: dup_1&dup_2 pair
	    block_1 = dup_1;
	    block_2 = dup_2;
	    donorRatio = p_donorRatio[2][4];
	  } else {			// WHC: wrong
	    cout << "something is wrong in void conversion()\n";
	    exit(0);
	  }

	  p = rand() / ((float) RAND_MAX + 1);
	  if (p < donorRatio) {
	    donor = block_1;
	    receptor = block_2;
	  } else {
	    donor = block_2;
	    receptor = block_1;
	  }
	  chrDonor = i;
	  chrReceptor = i;

	}
      } else if (chr1->b == 5) {
	// WHC: lost dup_1, but can still have IGC between ori and dup_2
	// WHC: forgot to add at the beginning

	cout << "Should use conversion_for_phaseVI()!\n";
	exit(0);
	int block_1 = ori_index;
	int block_2 = dup_2;
	float donorRatio = p_donorRatio[0][4];

	if (p < (2 * probability * BLOCKLENGTH * sameDifIGC)) {
	  IGC = 1;
	  // Determines which block will be the donor and which will be the receptor
	  p = rand() / ((float) RAND_MAX + 1);
	  if (p < donorRatio) {
	    donor = block_1;
	    receptor = block_2;
	  } else {
	    donor = block_2;
	    receptor = block_1;
	  }
	  chrDonor = i;
	  chrReceptor = i;
	}
      }
    }

    // WHC: not changed, remains the same; I assume it will work :)
    
    if(IGC == 1) {
      chr1 = pointer[pres][chrDonor];
      chr2 = pointer[pres][chrReceptor];
      // Determines conversion initiation point
      junction = (int) (rand() % (BLOCKLENGTH));
      // Determines gene conversion tract length through function tractpql()
      // WHC: what if this lies outside of a block???
      // WHC: it just doesn't copy the mutation outside this particular block; because mutation[BLOCK][mutation_index]
      tractlength = tractpql(meanTractLength);
      // If tractlength is an odd number, we must correct junction in order for gene conversion to be balanced
      if(tractlength%2 != 0){
	p = rand() / ((float) RAND_MAX + 1);
	if(p < 0.5){
	  junction += 1;
	}
      }
	
      // Tests for MEPS, 100% identity tract to allow conversion
      vald0 = location(junction - meps/2, pres, chrDonor, donor);
      valdf = location(junction + meps/2, pres, chrDonor, donor);
      valr0 = location(junction - meps/2, pres, chrReceptor, receptor);
      valrf = location(junction + meps/2, pres, chrReceptor, receptor);
      val = (valdf - vald0)-(valrf - valr0); // val is the number of positions to shift
      eqmut = 0;
      // If the sequences have the same number of mutations in the total identity region
      if (val == 0){
	for (k = 0; k < (valdf - vald0); k++) {
	  if (chr2[0].mutation[receptor][k + valr0] == chr1[0].mutation[donor][k + vald0]){
	    eqmut++;
	  }
	}
	// If the sequences have the same mutations in the total identity region
	if (eqmut == (valdf - vald0)){
	  // Determines location of the delimiting mutations for each block (donor and receptor)
	  vald0 = location(junction - 0.5 * tractlength, pres, chrDonor, donor);
	  valdf = location(junction + 0.5 * tractlength, pres, chrDonor, donor);
	  valr0 = location(junction - 0.5 * tractlength, pres, chrReceptor, receptor);
	  valrf = location(junction + 0.5 * tractlength, pres, chrReceptor, receptor);
	  val = (valdf - vald0)-(valrf - valr0); // val is the number of positions to shift
	  // Calculates differences between both sequences
	  differences = 0;
	  //Corrected bug; it used to say k<=
	  for (k = vald0; k < valdf; k++) {
	    otherk = location(chr1[0].mutation[donor][k], pres, chrReceptor, receptor);
	    if (chr2[0].mutation[receptor][otherk] != chr1[0].mutation[donor][k]){
	      differences++;
	    }
	  }
	  //Corrected bug; it used to say k<=
	  for (k = valr0; k < valrf; k++) {
	    otherk = location(chr2[0].mutation[receptor][k], pres, chrDonor, donor);
	    if (chr1[0].mutation[donor][otherk] != chr2[0].mutation[receptor][k]){
	      differences++;
	    }
	  }
	  // If differences do not exceed the similarity threshold for all the tractlength (MESH 2.0)
	  if(differences <= tractlength*(1-(similarityInConvTract/100))){
	    // Change receptor mutations in chr2[0].mutation[receptor][] array
	    if (val > 0) {
	      for (k = chr2[0].mpb[receptor] + val - 1; k >= valrf + val; k--) {
		chr2[0].mutation[receptor][k] = chr2[0].mutation[receptor][k - val];
	      }
	    }
	    if (val < 0) {
	      for (k = valrf + val; k < chr2[0].mpb[receptor] + val; k++) {
		chr2[0].mutation[receptor][k] = chr2[0].mutation[receptor][k - val];
	      }
	    }
	    for (k = 0; k < (valdf - vald0); k++) {
	      chr2[0].mutation[receptor][k + valr0] = chr1[0].mutation[donor][k + vald0];
	    }
	    // Change receptor chr2[0].mpb[] array
	    chr2[0].mpb[receptor] += val;
	  }
	}
      }
    }
  }
}


/* ================================================================ */
// WHC: conversion() for phaseVI()
void conversion_for_phaseVI (float probability, int t, int i, int pres, float (*p_donorRatio)[5], float sameDifIGC) {
  int k, partner, donor, receptor, chrDonor, chrReceptor, IGC, vald0, valdf, valr0, valrf, val, junction, tractlength, eqmut, otherk, differences;

  float p;
  chrom * chr1, * chr2;
  chr1 = pointer[pres][i];
  IGC = 0;



  
  // If the chr has duplicated block
  //  if (chr1[0].b == 3) {
  // WHC: as long as #blocks >= 3
  // WHC: REMEMBER, #blocks also determine which 2 pairs of duplications could be selected
  
  if (chr1->b >= 5) {
    if((sameDifIGC!=1) && (IGCmatrix[i][t]==true)){ // IGC on different chroms
      IGC = 1;

      if (i % 2 == 0) {
	partner = (i + 1);
      } else {
	partner = (i - 1);
      }			
      chr2 = pointer[pres][partner];

      /* ================================================================================================ */
      // beginning of new donor/receptor selection method
      
      // WHC: randomly choose a pair of duplicated blocks, depending on the number of blocks chrom has
      // WHC: also, return value of donorRatio
      int block_1, block_2;
      float donorRatio;
      
      if (chr1->b == 5 || chr2->b == 5) { // WHC: pick 2 pairs
	// WHC: REMINDER!! Probably should be if (chr1->b == 5 && chr2->b == 5) and add another if ( || )
	// WHC: IGC could only happend between ori and dup_2
	block_1 = ori_index;
	block_2 = dup_2;
	//	donorRatio = p_donorRatio[0][2]; WHC: mistake!!
	donorRatio = p_donorRatio[0][4];
      } else if (chr1->b == 6 && chr2->b == 6) {	// WHC: could pick ori, dup_1, dup_2 blocks; means already in phaseIV(), at least 3 blocks for any chrom
	// WHC: randomly generate 0, 1, 2
	// WHC: which correspond to ori&dup_1, ori&dup_2 and dup_1&dup_2 pairs

	int which_pair = pick_a_pair(1);
	//	int which_pair = rand() % 3;

	if (which_pair == 0) {	// WHC: ori&dup_1 pair
	  block_1 = ori_index;
	  block_2 = dup_1;
	  donorRatio = p_donorRatio[0][2];
	} else if (which_pair == 1) { // WHC: ori&dup_2 pair
	  block_1 = ori_index;
	  block_2 = dup_2;
	  donorRatio = p_donorRatio[0][4];
	} else if (which_pair == 2) { // WHC: dup_1&dup_2 pair
	  block_1 = dup_1;
	  block_2 = dup_2;
	  donorRatio = p_donorRatio[2][4];
	} else {			// WHC: wrong
	  cout << "something is wrong in void conversion()\n";
	  exit(0);
	}
      } else {
	cout << "wrong on void conversion_for_phaseVI() 2\n";
	exit(0);
      }

      // Determines which block will be the donor and which will be the receptor
      p = rand() / ((float) RAND_MAX + 1);
      if (p < donorRatio) {
	donor = block_1;
	receptor = block_2;
      } else {
	donor = block_2;
	receptor = block_1;
      }

      //      if (chr2->b == 5 || chr1->b == 5) {; not necessary anymore
	p = rand() / ((float) RAND_MAX + 1); // WHC: randomly choose from which chromosome to which
	if (p < 0.5) {
	  chrDonor = i;
	  chrReceptor = partner;
	} else {
	  chrDonor = partner;
	  chrReceptor = i;
	}


      // end of new donor/receptor selection method
      /* ================================================================================================ */
      
    } else {
      // WHC: if IGCmatrix[i][t] == false, this means that no IGC between chroms will happen, due to void genealogy()
      // WHC: but IGC could still happen within chrom
      p = rand() / ((float) RAND_MAX + 1);
      // WHC: should distinguish between 2 duplications and 3 duplications
      
      if (chr1->b == 5) {
	int block_1 = ori_index;
	int block_2 = dup_2;
	float donorRatio = p_donorRatio[0][4];

	if (p < (2 * probability * BLOCKLENGTH * sameDifIGC)) {
	  IGC = 1;
	  // Determines which block will be the donor and which will be the receptor
	  p = rand() / ((float) RAND_MAX + 1);
	  if (p < donorRatio) {
	    donor = block_1;
	    receptor = block_2;
	  } else {
	    donor = block_2;
	    receptor = block_1;
	  }
	  chrDonor = i;
	  chrReceptor = i;
	}
      } else if (chr1->b == 6) {
	if (p < (3 * probability * BLOCKLENGTH * sameDifIGC)) { // WHC: should this be 3? Because we have 3 blocks...
	  IGC = 1;
	  // Determines which block will be the donor and which will be the receptor

	  int which_pair = pick_a_pair(1);
	  //	  int which_pair = rand() % 3;
	  int block_1, block_2;
	  float donorRatio;
	  if (which_pair == 0) {	// WHC: ori&dup_1 pair
	    block_1 = ori_index;
	    block_2 = dup_1;
	    donorRatio = p_donorRatio[0][2];
	  } else if (which_pair == 1) { // WHC: ori&dup_2 pair
	    block_1 = ori_index;
	    block_2 = dup_2;
	    donorRatio = p_donorRatio[0][4];
	  } else if (which_pair == 2) { // WHC: dup_1&dup_2 pair
	    block_1 = dup_1;
	    block_2 = dup_2;
	    donorRatio = p_donorRatio[2][4];
	  } else {			// WHC: wrong
	    cout << "something is wrong in void conversion()\n";
	    exit(0);
	  }

	  p = rand() / ((float) RAND_MAX + 1);
	  if (p < donorRatio) {
	    donor = block_1;
	    receptor = block_2;
	  } else {
	    donor = block_2;
	    receptor = block_1;
	  }
	  chrDonor = i;
	  chrReceptor = i;

	}
      }
    }

    // WHC: not changed, remains the same; I assume it will work :)
    
    if(IGC == 1) {
      //      cout << "IGC should not equal to 1\n";
      //      exit(0);
      chr1 = pointer[pres][chrDonor];
      chr2 = pointer[pres][chrReceptor];
      // Determines conversion initiation point
      junction = (int) (rand() % (BLOCKLENGTH));
      // Determines gene conversion tract length through function tractpql()
      // WHC: what if this lies outside of a block???
      // WHC: it just doesn't copy the mutation outside this particular block; because mutation[BLOCK][mutation_index]
      tractlength = tractpql(meanTractLength);
      // If tractlength is an odd number, we must correct junction in order for gene conversion to be balanced
      if(tractlength%2 != 0){
	p = rand() / ((float) RAND_MAX + 1);
	if(p < 0.5){
	  junction += 1;
	}
      }
	
      // Tests for MEPS, 100% identity tract to allow conversion
      vald0 = location(junction - meps/2, pres, chrDonor, donor);
      valdf = location(junction + meps/2, pres, chrDonor, donor);
      valr0 = location(junction - meps/2, pres, chrReceptor, receptor);
      valrf = location(junction + meps/2, pres, chrReceptor, receptor);
      val = (valdf - vald0)-(valrf - valr0); // val is the number of positions to shift
      eqmut = 0;
      // If the sequences have the same number of mutations in the total identity region
      if (val == 0){
	for (k = 0; k < (valdf - vald0); k++) {
	  if (chr2[0].mutation[receptor][k + valr0] == chr1[0].mutation[donor][k + vald0]){
	    eqmut++;
	  }
	}
	// If the sequences have the same mutations in the total identity region
	if (eqmut == (valdf - vald0)){
	  // Determines location of the delimiting mutations for each block (donor and receptor)
	  vald0 = location(junction - 0.5 * tractlength, pres, chrDonor, donor);
	  valdf = location(junction + 0.5 * tractlength, pres, chrDonor, donor);
	  valr0 = location(junction - 0.5 * tractlength, pres, chrReceptor, receptor);
	  valrf = location(junction + 0.5 * tractlength, pres, chrReceptor, receptor);
	  val = (valdf - vald0)-(valrf - valr0); // val is the number of positions to shift
	  // Calculates differences between both sequences
	  differences = 0;
	  //Corrected bug; it used to say k<=
	  for (k = vald0; k < valdf; k++) {
	    otherk = location(chr1[0].mutation[donor][k], pres, chrReceptor, receptor);
	    if (chr2[0].mutation[receptor][otherk] != chr1[0].mutation[donor][k]){
	      differences++;
	    }
	  }
	  //Corrected bug; it used to say k<=
	  for (k = valr0; k < valrf; k++) {
	    otherk = location(chr2[0].mutation[receptor][k], pres, chrDonor, donor);
	    if (chr1[0].mutation[donor][otherk] != chr2[0].mutation[receptor][k]){
	      differences++;
	    }
	  }
	  // If differences do not exceed the similarity threshold for all the tractlength (MESH 2.0)
	  if(differences <= tractlength*(1-(similarityInConvTract/100))){
	    // Change receptor mutations in chr2[0].mutation[receptor][] array
	    if (val > 0) {
	      for (k = chr2[0].mpb[receptor] + val - 1; k >= valrf + val; k--) {
		chr2[0].mutation[receptor][k] = chr2[0].mutation[receptor][k - val];
	      }
	    }
	    if (val < 0) {
	      for (k = valrf + val; k < chr2[0].mpb[receptor] + val; k++) {
		chr2[0].mutation[receptor][k] = chr2[0].mutation[receptor][k - val];
	      }
	    }
	    for (k = 0; k < (valdf - vald0); k++) {
	      chr2[0].mutation[receptor][k + valr0] = chr1[0].mutation[donor][k + vald0];
	    }
	    // Change receptor chr2[0].mpb[] array
	    chr2[0].mpb[receptor] += val;
	  }
	}
      }
    }
  } else { cout << "mistake happens.\n"; exit(0); }
}

/* ================================================================ */

//////////////////////
////  STATISTICS  ////
//////////////////////

void statistics(int prev, bool does_print) {
  int j, o;
  float * resultsSample;
  //  float results_pairwise_divergence_between;
  
  // WHC: FSL() will delete fixed/lost mutations
  // WHC: and make a new muttable[]

  /* TESTING PURPOSE */
  /*
  cout << "The MutCount is " << MutCount << '\n';
  int num_of_true_in_multihit = 0, num_of_false_in_multihit = 0;
  for (int i = 0; i < BLOCKLENGTH; ++i) {
    if (multihit[i] == true) {
      num_of_true_in_multihit++;
    } else if (multihit[i] == false) {
      num_of_false_in_multihit++;
    } else {
      cout << "WRONG WRONG\n";
      exit(0);
    }
  }
  cout << "The number of true in multihit is " << num_of_true_in_multihit << " and false is " << num_of_false_in_multihit << endl;

  int num_of_true_in_multihit_single = 0, num_of_false_in_multihit_single = 0;
  for (int i = 0; i < BLOCKLENGTH; ++i) {
    if (multihit_single[i] == true) {
      num_of_true_in_multihit_single++;
    } else if (multihit_single[i] == false) {
      num_of_false_in_multihit_single++;
    } else {
      cout << "WRONG WRONG\n";
      exit(0);
    }
  }
  cout << "The number of true in multihit_single is " << num_of_true_in_multihit_single << " and false is " << num_of_false_in_multihit_single << endl;
  cout << "Total number of true is " << num_of_true_in_multihit + num_of_true_in_multihit_single << endl;
  */
  
  if (duFreq_2 == false) {
    FSL(prev); // Creating the summary vector fixedLostForAll
    // INFORMATION RECOVERED FROM SAMPLES
  } else {
    FSL_since_phaseIV(prev);
  }
  
  for (o = 0; o < numofsamples; o++) {
    SamplingIndividuals(sampleN[o]);

    // CALCULATES THE SFS FOR THE SAMPLE BLOCK BY BLOCK
    for (j = 0; j < B; j++) {
      resultsSample = SiteFrequencySpectrumPrint(prev, j, sampleN[o], does_print);
      samplefile[j][o][0] << resultsSample[1] << " "; // the results array keeps (S,pi) but for the sake of tradition
      samplefile[j][o][1] << resultsSample[0] << " "; // we save it as (pi,S)
    }

    /* =============================================================== */
    // It works!!! Their values are consistently higher than mine, as they used integer divisons...

    /*
    // WHC: just to test the differences between these 2 methods of calculating pairwise divergence
    resultsSample = SiteFrequencySpectrumPrint(prev, 0, sampleN[o], does_print);
    cout << "from their calculation: " << resultsSample[1] << '\n';

    // WHC: my calculations
    
    float temp_results = 0.0;
    for (int c_1 = 0; c_1 < 2 * sampleN[o]; ++c_1) {
      for (int c_2 = 0; c_2 < 2 * sampleN[o]; ++c_2) {
	if (c_1 == c_2) {
	  if (number_of_nucleotide_diff(0, 0, sample[c_1], sample[c_2], prev) != 0) { // another test
	    cout << "hi, something wrong with the number_of_nucleotide_diff() function!\n";
	    exit(0);
	  }
	  continue;
	}
	temp_results += number_of_nucleotide_diff(0, 0, sample[c_1], sample[c_2], prev);
      }
    }
    // mine results and theirs always have a difference by factor 2n/(2n-1)
    cout << "from my calculation: " << temp_results / (2 * sampleN[o] * 2 * sampleN[o]) << '\n';
    cout << "from my adjusted calculation: " << temp_results / (2 * sampleN[o] * (2 * sampleN[o] - 1)) << '\n';
    
    */
    /* ================================================================ */
    
    // WHC: for calculating the pairwise-divergence between blockA and blockB
    // WHC: will just try to calculate blockA = 0, blockB = 3 for now
    if (duFreq == false) {
      // do nothing
      // NO, should cout << 0!
      for (int k = 0; k < 3; ++k) {
	interblock_divergence[k] << '0' << " ";
      }
    } else if (duFreq == true && duFreq_2 == false) {
      cout << "between original and dup_1: " << pairwise_divergence_between(prev, sampleN[o], 0, 3) << '\n';
      
      // between original and dup_1
      interblock_divergence[0] << pairwise_divergence_between(prev, sampleN[o], 0, 3) << " ";
      // between original and dup_2
      interblock_divergence[1] << 0 << " ";
      // between dup_1 and dup_2
      interblock_divergence[2] << 0 << " ";
      
    } else if (duFreq_2 == true && loseFreq == false) {
      cout << "between original and dup_1: " << pairwise_divergence_between(prev, sampleN[o], 0, 3) << '\n';
      cout << "between original and dup_2: " << pairwise_divergence_between(prev, sampleN[o], 0, 5) << '\n';
      cout << "between dup_1 and dup_2: " << pairwise_divergence_between(prev, sampleN[o], 3, 5) << '\n';

      // between original and dup_1
      interblock_divergence[0] << pairwise_divergence_between(prev, sampleN[o], 0, 3) << " ";
      // between original and dup_2
      interblock_divergence[1] << pairwise_divergence_between(prev, sampleN[o], 0, 5) << " ";
      // between dup_1 and dup_2
      interblock_divergence[2] << pairwise_divergence_between(prev, sampleN[o], 3, 5) << " ";
      
    } else if (loseFreq == true) {
      cout << "should not be here.\n";
      exit(0);
    } else {
      cout << "should not be here.\n";
      exit(0);
    }

  }
  // CALCULATES THE NUMBER OF PRIVATE & SHARED MUTATIONES BETWEEN BLOCKS 0 & 2
  DivergenceForAll(0, 3, prev);
  DivergenceForAll(0, 5, prev);
  DivergenceForAll(3, 5, prev);
}

/* ================================================================ */
// WHC: new statistics() for phaseVI; need to creat SiteFrequencySpectrumPrint_for_phaseVI()

void statistics_for_phaseVI(int prev, bool does_print) {
  int j, o;
  float * resultsSample;

  // WHC: FSL() will delete fixed/lost mutations
  // WHC: and make a new muttable[]

  // WHC: as location() and muFrequencyIntWholePopAndSample() should work with block = 2, I think FSL() will work too
  //  FSL(prev); // Creating the summary vector fixedLostForAll
  // INFORMATION RECOVERED FROM SAMPLES

  /* TESTING PURPOSE */
  /*
  cout << "The MutCount is " << MutCount << '\n';
  int num_of_true_in_multihit = 0, num_of_false_in_multihit = 0;
  for (int i = 0; i < BLOCKLENGTH; ++i) {
    if (multihit[i] == true) {
      num_of_true_in_multihit++;
    } else if (multihit[i] == false) {
      num_of_false_in_multihit++;
    } else {
      cout << "WRONG WRONG\n";
      exit(0);
    }
  }
  cout << "The number of true in multihit is " << num_of_true_in_multihit << " and false is " << num_of_false_in_multihit << endl;

  int num_of_true_in_multihit_single = 0, num_of_false_in_multihit_single = 0;
  for (int i = 0; i < BLOCKLENGTH; ++i) {
    if (multihit_single[i] == true) {
      num_of_true_in_multihit_single++;
    } else if (multihit_single[i] == false) {
      num_of_false_in_multihit_single++;
    } else {
      cout << "WRONG WRONG\n";
      exit(0);
    }
  }
  cout << "The number of true in multihit_single is " << num_of_true_in_multihit_single << " and false is " << num_of_false_in_multihit_single << endl;
  cout << "Total number of true is " << num_of_true_in_multihit + num_of_true_in_multihit_single << endl;

  */

  
  FSL_since_phaseIV(prev);
  
  for (o = 0; o < numofsamples; o++) {
    SamplingIndividuals(sampleN[o]);

    // CALCULATES THE SFS FOR THE SAMPLE BLOCK BY BLOCK
    for (j = 0; j < B; j++) {
      resultsSample = SiteFrequencySpectrumPrint_for_phaseVI(prev, j, sampleN[o], does_print);
      samplefile[j][o][0] << resultsSample[1] << " "; // the results array keeps (S,pi) but for the sake of tradition
      samplefile[j][o][1] << resultsSample[0] << " "; // we save it as (pi,S)
    }


    
    // WHC: for calculating the pairwise-divergence between blockA and blockB
    // WHC: will just try to calculate blockA = 0, blockB = 3 for now
    if (loseFreq == false) { cout << "probably should not use this.\n"; exit(0); }
    
    cout << "between original and dup_1: " << pairwise_divergence_between(prev, sampleN[o], 0, 3) << '\n';
    interblock_divergence[0] << pairwise_divergence_between(prev, sampleN[o], 0, 3) << " ";
    
    cout << "between original and dup_2: " << pairwise_divergence_between(prev, sampleN[o], 0, 5) << '\n';
    interblock_divergence[1] << pairwise_divergence_between(prev, sampleN[o], 0, 5) << " ";
    
    cout << "between dup_1 and dup_2: " << pairwise_divergence_between(prev, sampleN[o], 5, 3) << '\n';
    interblock_divergence[2] << pairwise_divergence_between(prev, sampleN[o], 5, 3) << " ";
  }
  // CALCULATES THE NUMBER OF PRIVATE & SHARED MUTATIONES BETWEEN BLOCKS 0 & 2
  DivergenceForAll(0, 3, prev);
  DivergenceForAll(0, 5, prev);
  DivergenceForAll(3, 5, prev);
}

/* ================================================================ */

/////////////////////////////////////////////////////////////////////
////  FSL (FIXED-SEGREGATING-LOST): EXTRACTS INFO FROM MUTTABLE  ////
/////////////////////////////////////////////////////////////////////

// WHC: FSL() is changed to accomodate phaseIV()
void FSL(int hh) {
  int m, otherm=0, mm;

  for (m = 0; m < MutCount; m++) {
    muttable[m].frequency = muFrequencyIntWholePopAndSample(hh, muttable[m].position, muttable[m].block, N);
    // WHC: What if this block is lost, should I use (N * 2) or the #chroms carrying this block????
    // WHC: the previous question may be vitally important, since this could cause a continuously increase of multihit[]==true

    // WHC: IMPORTANT need to doulbe check here!!!!???
    
    // muttable[m].frequency = ((float) muttable[m].frequency) / (N * 2);
    if (muttable[m].block == 3 && loseFreq == false) {
      muttable[m].frequency = ((float) muttable[m].frequency) / (DupliFreq(hh, 3, N));
    } else if (muttable[m].block == 3 && loseFreq == true) {
      muttable[m].frequency = ((float) muttable[m].frequency) / (DupliFreq_for_phaseVI(hh, 3, N));
    } else if (muttable[m].block == 5 && loseFreq == false) {
      // in phaseIV (when block_4 not fixed) and in phaseV (which should always be 2*N
      muttable[m].frequency = ((float) muttable[m].frequency) / (DupliFreq(hh, 5, N));
    } else if (muttable[m].block == 5 && loseFreq == true) {
      muttable[m].frequency = ((float) muttable[m].frequency) / (DupliFreq(hh, 5, N));
      if (2 * N != DupliFreq(hh, 5, N)) {
	cout << "This 2 values should equal!\n";
	exit(0);
      }
    } else if (muttable[m].block == 0 || muttable[m].block == 1) {
      muttable[m].frequency = ((float) muttable[m].frequency) / (N * 2);
    } else if (muttable[m].block == 2 || muttable[m].block == 4) {
      // do nothing
      continue;
    } else {
      cout << "INcorrect\n";
      exit(0);
    }
    
  }
  for (m = 0, mm = 0; m < MutCount; m++) {
    // SET MULTIHIT TO FALSE FOR MUTATIONS THAT HAVE BEEN LOST

    // WHC: this -9 labels that this has been added to tempMutCount[], so that it won't be added twice
    if (muttable[m].block == -9 || muttable[m].block == 2 || muttable[m].block == 4) {
      continue;
    }

    if (muttable[m].frequency == 0 && muttable[m].block == 1) {
      // WHC: this should not be correct???? single-block should have its own multihit vector????
      
      multihit_single[muttable[m].position] = false;
    }
    if (muttable[m].frequency == 0 && muttable[m].block == 0 && duFreq == false) {
      multihit[muttable[m].position] = false;
    }
    // IF MUTATION IS IN THE SINGLE-COPY BLOCK
    if (muttable[m].block == 1) {
      // SEGREGATING
      if (muttable[m].frequency > 0 && muttable[m].frequency < 1) {
	temporalmuttable[mm].block = muttable[m].block;
	temporalmuttable[mm].frequency = muttable[m].frequency;
	temporalmuttable[mm].position = muttable[m].position;
	mm++;

      } else if (muttable[m].frequency == 1) {      // FIXED
	EraseFixedMutations(muttable[m].position, muttable[m].block, hh);
      }

    }// IF THE MUTATION IS IN EITHER BLOCK 0 0R 2
     // WHC: OR BLOCK 4
    // WHC: in block 0, 3, 5
    else {
      // IF THE DUPLICATION HAS NOT YET OCCURRED, TREATS BLOCK 0 AS SINGLE-COPY
      // WHC: if duFreq = false, then duFreq_2 = false must be true
      if (duFreq == false) {
	// SEGREGATING
	if (muttable[m].frequency > 0 && muttable[m].frequency < 1) {
	  temporalmuttable[mm].block = muttable[m].block;
	  temporalmuttable[mm].frequency = muttable[m].frequency;
	  temporalmuttable[mm].position = muttable[m].position;
	  mm++;
	}
	// FIXED
	if (muttable[m].frequency == 1) {
	  EraseFixedMutations(muttable[m].position, muttable[m].block, hh);
	}
      }  else if ((duFreq == true) && (duFreq_2 == false)) {
        // IF THE DUPLICATION HAS OCCURRED
	// WHC: but before phaseIV()
	
	// CHECKS THE BLOCK IN WHICH IT HAS OCCURRED
	// AND LOOKS FOR THE POSITION OF THE MUTATION IN MUTTABLE FOR THE OTHER BLOCK
	if (muttable[m].block == 3) {
	  otherm = SearchMutation(0, muttable[m].position, MutCount);
	}
	if (muttable[m].block == 0) {
	  otherm = SearchMutation(3, muttable[m].position, MutCount);
	}
	// SEGREGATING
	if (muttable[m].frequency > 0 && muttable[m].frequency < 1) {
	  temporalmuttable[mm].block = muttable[m].block;
	  temporalmuttable[mm].frequency = muttable[m].frequency;
	  temporalmuttable[mm].position = muttable[m].position;
	  mm++;
	  // IN CASE THE MUTATION IS NOT PRESENT IN THE OTHER BLOCK, IT IS COPIED TO MUTTABLE ANYWAY
	  if (muttable[otherm].frequency == 0) {
	    temporalmuttable[mm].block = muttable[otherm].block;
	    temporalmuttable[mm].frequency = muttable[otherm].frequency;
	    temporalmuttable[mm].position = muttable[otherm].position;
	    mm++;
	  }
	}
	// FIXED
	if (muttable[m].frequency == 1) {
	  // IF THE DUPLICATION IS FIXED IN THE OTHER BLOCK
	  if (muttable[otherm].frequency == 1) {
	    // WHC: why don't erase mutaitons in muttable[otherm]???
	    // WHC: because otherm will go through this later?
	    // WHC: and when going through mutation[otherm], because its partner is erased, won't enter this block(?)
	    // WHC: I think it should erase mutation[otherm] too!
	    // WHC: NO! Do not need to; because mutation[m].frequency will NOT change!
	    EraseFixedMutations(muttable[m].position, muttable[m].block, hh);
	  }// IF THE DUPLICATION IS SEGREGATING OR NOT PRESENT IN THE OTHER BLOCK
	  else {
	    temporalmuttable[mm].block = muttable[m].block;
	    temporalmuttable[mm].frequency = muttable[m].frequency;
	    temporalmuttable[mm].position = muttable[m].position;
	    mm++;
	    // IN CASE THE MUTATION IS NOT PRESENT IN THE OTHER BLOCK, IT IS COPIED TO MUTTABLE ANYWAY
	    if (muttable[otherm].frequency == 0) {
	      // WHC: because there is NO if (muttable[m].frequency == 0) {} clause, this copy here won't cause copying twice
	      temporalmuttable[mm].block = muttable[otherm].block;
	      temporalmuttable[mm].frequency = muttable[otherm].frequency;
	      temporalmuttable[mm].position = muttable[otherm].position;
	      mm++;
	    }
	  }
	}
	// LOST (TO CONSIDER IT LOST IT WOULD HAVE TO BE LOST IN BOTH BLOCKS)
	if (muttable[m].frequency == 0) {
	  if (muttable[otherm].frequency == 0) {
	    multihit[muttable[m].position] = false;
	  }
	}
      } else if (duFreq_2 == true) {
	cout << "Should use FSL_since_phaseIV().\n";
	exit(0);
	int otherm_1 = 0, otherm_2 = 0;

	/* ================================================================ */
	/* WHC: this clause is a modificate of if (duFreq == true && duFreq_2 == false) */
	
        // IF THE SECONDE DUPLICATION HAS OCCURRED
	// WHC: within or after phaseIV() (probabily will write one specific for phaseVI()...
	
	// CHECKS THE BLOCK IN WHICH IT HAS OCCURRED
	// AND LOOKS FOR THE POSITION OF THE MUTATION IN MUTTABLE FOR THE OTHER BLOCK
	if (muttable[m].block == 2) {
	  otherm_1 = SearchMutation(0, muttable[m].position, MutCount);
	  otherm_2 = SearchMutation(4, muttable[m].position, MutCount);
	} else if (muttable[m].block == 0) {
	  otherm_1 = SearchMutation(2, muttable[m].position, MutCount);
	  otherm_2 = SearchMutation(4, muttable[m].position, MutCount);
	} else if (muttable[m].block == 4) {
	  otherm_1 = SearchMutation(0, muttable[m].position, MutCount);
	  otherm_2 = SearchMutation(2, muttable[m].position, MutCount);
	} else if (muttable[m].block == 3) {
	  // WHC: we just don't care about what's going on in 3
	  cout << "Should already skipped block == 3\n";
	  exit(0);
	} else {
	  cout << "hi, something wrong in FSL() when considering phaseIV().\n";
	  exit(0);
	}

	// SEGREGATING
	
	// WHC: consider a senario with frequencies: 1, 0.5, 0 for m, otherm_1, otherm_2 !!!
	
	if (muttable[m].frequency > 0 && muttable[m].frequency < 1) {
	  temporalmuttable[mm].block = muttable[m].block;
	  temporalmuttable[mm].frequency = muttable[m].frequency;
	  temporalmuttable[mm].position = muttable[m].position;
	  mm++;

	  //	  if (muttable[otherm_1].frequency == 0 || muttable[otherm_2].frequency == 0) {
	  // WHC: consider a senario with frequencies: 1, 0.5, 0 for m, otherm_1, otherm_2 !!!
	  if (muttable[otherm_1].frequency == 0 && muttable[otherm_1].block != -9) {
	    temporalmuttable[mm].block = muttable[otherm_1].block;
	    temporalmuttable[mm].frequency = muttable[otherm_1].frequency;
	    temporalmuttable[mm].position = muttable[otherm_1].position;
	    mm++;

	    muttable[otherm_1].block = -9;
	  }

	  if (muttable[otherm_2].frequency == 0 && muttable[otherm_2].block != -9) {
	    temporalmuttable[mm].block = muttable[otherm_2].block;
	    temporalmuttable[mm].frequency = muttable[otherm_2].frequency;
	    temporalmuttable[mm].position = muttable[otherm_2].position;
	    mm++;

	    // WHC: a label, so that it will skip at the beginning

	    muttable[otherm_2].block = -9;
	  }

	  /* WHC: will merge them together, because...
	  // IN CASE THE MUTATION IS NOT PRESENT IN THE OTHER BLOCK, IT IS COPIED TO MUTTABLE ANYWAY
	  if (muttable[otherm_1].frequency == 0) {
	    temporalmuttable[mm].block = muttable[otherm_1].block;
	    temporalmuttable[mm].frequency = muttable[otherm_1].frequency;
	    temporalmuttable[mm].position = muttable[otherm_1].position;
	    mm++;

	    // WHC: a label, so that it will skip at the beginning
	    muttable[otherm_1].block = -9;
	  }

	  // WHC: can NOT do this! because it will copy them more than once!!!
	  // WHC: if 0 < muttable[m] < 1, and 0 < muttable[otherm_1] < 1, and muttable[otherm_2] == 0
	  // WHC: then, muttable[otherm_2] will be added to temporalmuttable[] TWICE!!!

	  if (muttable[otherm_2].frequency == 0) {
	    temporalmuttable[mm].block = muttable[otherm_2].block;
	    temporalmuttable[mm].frequency = muttable[otherm_2].frequency;
	    temporalmuttable[mm].position = muttable[otherm_2].position;
	    mm++;

	    // WHC: a label, so that it will skip at the beginning
	    muttable[otherm_2].block = -9;

	  }
	  */

	  // WHC: maybe, I could lable the original muttable[otherm_1].block = -9, and skip it using continue clause??

	} else if (muttable[m].frequency == 1) {	// FIXED
	  // IF THE DUPLICATION IS FIXED IN THE OTHER BLOCK
	  if ((muttable[otherm_1].frequency == 1) && (muttable[otherm_2].frequency == 1)) {
	    // WHC: all 3 are fixed
	    
	    // WHC: why don't erase mutaitons in muttable[otherm]???
	    // WHC: because otherm will go through this later?
	    // WHC: and when going through mutation[otherm], because its partner is erased, won't enter this block(?)
	    // WHC: I think it should erase mutation[otherm] too!
	    // WHC: NO! Do not need to; because mutation[m].frequency will NOT change!
	    
	    EraseFixedMutations(muttable[m].position, muttable[m].block, hh); // WHC: only erase muttable[m], not muttable[otherm]

	    // WHC: should I change multihit[] = false???
	    // WHC: NO. EraseFixedMutations() already did that for you
	    
	  }// IF THE DUPLICATION IS SEGREGATING OR NOT PRESENT IN THE OTHER BLOCK
	  else {
	    temporalmuttable[mm].block = muttable[m].block;
	    temporalmuttable[mm].frequency = muttable[m].frequency;
	    temporalmuttable[mm].position = muttable[m].position;
	    mm++;
	    // IN CASE THE MUTATION IS NOT PRESENT IN THE OTHER BLOCK, IT IS COPIED TO MUTTABLE ANYWAY
	    // WHC consider frequency 0.5, 1, 0, so should be && but not || (because 0.5 will be added twice!)

	    /*

	    if (muttable[otherm_1].frequency == 0 && muttable[otherm_2].frequency == 0) {
	      // WHC: because there is NO if (muttable[m].frequency == 0) {} clause, this copy here won't cause copying twice
	      temporalmuttable[mm].block = muttable[otherm_1].block;
	      temporalmuttable[mm].frequency = muttable[otherm_1].frequency;
	      temporalmuttable[mm].position = muttable[otherm_1].position;
	      mm++;

	      temporalmuttable[mm].block = muttable[otherm_2].block;
	      temporalmuttable[mm].frequency = muttable[otherm_2].frequency;
	      temporalmuttable[mm].position = muttable[otherm_2].position;
	      mm++;

	      // WHC: a label, so that it will skip at the beginning
	      muttable[otherm_1].block = -9;
	      muttable[otherm_2].block = -9;

	    }

	    */
	    if (muttable[otherm_1].frequency == 0 && muttable[otherm_1].block != -9) {
	    temporalmuttable[mm].block = muttable[otherm_1].block;
	    temporalmuttable[mm].frequency = muttable[otherm_1].frequency;
	    temporalmuttable[mm].position = muttable[otherm_1].position;
	    mm++;

	    muttable[otherm_1].block = -9;
	  }

	  if (muttable[otherm_2].frequency == 0 && muttable[otherm_2].block != -9) {
	    temporalmuttable[mm].block = muttable[otherm_2].block;
	    temporalmuttable[mm].frequency = muttable[otherm_2].frequency;
	    temporalmuttable[mm].position = muttable[otherm_2].position;
	    mm++;

	    // WHC: a label, so that it will skip at the beginning

	    muttable[otherm_2].block = -9;
	  }

	  
	  }
	} else if (muttable[m].frequency == 0) {	// LOST (TO CONSIDER IT LOST IT WOULD HAVE TO BE LOST IN BOTH BLOCKS)
	  if (muttable[otherm_1].frequency == 0 && muttable[otherm_2].frequency == 0) {
	    multihit[muttable[m].position] = false;
	  }
	} else {
	  cout << "something is wrong in FSL() here.\n";
	  exit(0);
	}

	
	/* WHC: this clause is a modificate of if (duFreq == true && duFreq_2 == false) */
	/* ================================================================ */
      } else {
	cout << "something is wrong in FSL\n";
	exit(0);
      }
    }
  }
  MutCount = mm;
  for (m = 0; m < MutCount; m++) {
    muttable[m].block = temporalmuttable[m].block;
    muttable[m].frequency = temporalmuttable[m].frequency;
    muttable[m].position = temporalmuttable[m].position;
  }
}


/* ======================================================================================================================= */

void FSL_since_phaseIV(int hh) {


  int m, otherm=0, mm;

  int make_false_count = 0;
  
  for (m = 0; m < MutCount; m++) {
    muttable[m].frequency = muFrequencyIntWholePopAndSample(hh, muttable[m].position, muttable[m].block, N);
    // WHC: What if this block is lost, should I use (N * 2) or the #chroms carrying this block????
    // WHC: the previous question may be vitally important, since this could cause a continuously increase of multihit[]==true

    // WHC: IMPORTANT need to doulbe check here!!!!???
    
    // muttable[m].frequency = ((float) muttable[m].frequency) / (N * 2);
    if (muttable[m].block == 3 && loseFreq == false) {
      muttable[m].frequency = ((float) muttable[m].frequency) / (DupliFreq(hh, 3, N));
    } else if (muttable[m].block == 3 && loseFreq == true) {
      muttable[m].frequency = ((float) muttable[m].frequency) / (DupliFreq_for_phaseVI(hh, 3, N));
    } else if (muttable[m].block == 5 && loseFreq == false) {
      // in phaseIV (when block_4 not fixed) and in phaseV (which should always be 2*N
      muttable[m].frequency = ((float) muttable[m].frequency) / (DupliFreq(hh, 5, N));
    } else if (muttable[m].block == 5 && loseFreq == true) {
      muttable[m].frequency = ((float) muttable[m].frequency) / (DupliFreq(hh, 5, N));
      if (2 * N != DupliFreq(hh, 5, N)) {
	cout << "This 2 values should equal!\n";
	exit(0);
      }
    } else if (muttable[m].block == 0 || muttable[m].block == 1) {
      muttable[m].frequency = ((float) muttable[m].frequency) / (N * 2);
    } else if (muttable[m].block == 2 || muttable[m].block == 4) {
      // do nothing
      // WHC: skip block 2 and 4, the two spacers
      continue;
    } else {
      cout << "INcorrect\n";
      exit(0);
    }

    if (muttable[m].frequency > 1) {
      cout << "The frequency cannot be over 1!\n";
      exit(0);
    }
    
  }
  for (m = 0, mm = 0; m < MutCount; m++) {
    // SET MULTIHIT TO FALSE FOR MUTATIONS THAT HAVE BEEN LOST

    // WHC: this -9 labels that this has been added to tempMutCount[], so that it won't be added twice
    if (muttable[m].block == 2 || muttable[m].block == 4 || muttable[m].block == -9) {
      continue;
    }

    // WHC: since alreay after phaseIV
    if (duFreq_2 != true) {
      cout << "WRONG! You shouldn't use this.\n";
      exit(0);
    }
    
    // IF MUTATION IS IN THE SINGLE-COPY BLOCK
    if (muttable[m].block == 1) {
      // SEGREGATING
      if (muttable[m].frequency > 0 && muttable[m].frequency < 1) {
	temporalmuttable[mm].block = muttable[m].block;
	temporalmuttable[mm].frequency = muttable[m].frequency;
	temporalmuttable[mm].position = muttable[m].position;
	mm++;

      } else if (muttable[m].frequency == 1) {
	EraseFixedMutations(muttable[m].position, muttable[m].block, hh);
      } else if (muttable[m].frequency == 0) {
	multihit_single[muttable[m].position] = false;
      }

    }// IF THE MUTATION IS IN EITHER BLOCK 0 0R 2
     // WHC: OR BLOCK 4
    // WHC: block 0, 3, 5
    else {
      // IF THE DUPLICATION HAS NOT YET OCCURRED, TREATS BLOCK 0 AS SINGLE-COPY
      // WHC: if duFreq = false, then duFreq_2 = false must be true

      if (duFreq_2 == true) {
	int otherm_1 = 0, otherm_2 = 0;

	/* ================================================================ */
	/* WHC: this clause is a modificate of if (duFreq == true && duFreq_2 == false) */
	
        // IF THE SECONDE DUPLICATION HAS OCCURRED
	// WHC: within or after phaseIV() (probabily will write one specific for phaseVI()...
	
	// CHECKS THE BLOCK IN WHICH IT HAS OCCURRED
	// AND LOOKS FOR THE POSITION OF THE MUTATION IN MUTTABLE FOR THE OTHER BLOCK
	if (muttable[m].block == 3) {
	  otherm_1 = SearchMutation(0, muttable[m].position, MutCount);
	  otherm_2 = SearchMutation(5, muttable[m].position, MutCount);
	} else if (muttable[m].block == 0) {
	  otherm_1 = SearchMutation(3, muttable[m].position, MutCount);
	  otherm_2 = SearchMutation(5, muttable[m].position, MutCount);
	} else if (muttable[m].block == 5) {
	  otherm_1 = SearchMutation(0, muttable[m].position, MutCount);
	  otherm_2 = SearchMutation(3, muttable[m].position, MutCount);
	} else if (muttable[m].block == 2 || muttable[m].block == 4) {
	  // WHC: we just don't care about what's going on in 3
	  cout << "Should already skipped block == 2 or 4\n";
	  exit(0);
	} else {
	  cout << "hi, something wrong in FSL() when considering phaseIV().\n";
	  exit(0);
	}

	if (muttable[m].position != muttable[otherm_1].position || muttable[m].position != muttable[otherm_2].position) {
	  cout << "check out here\n";
	  exit(0);
	}
	
	// WHC: fixed or lost
	if (muttable[m].frequency == 1 && muttable[otherm_1].frequency == 1 && muttable[otherm_2].frequency == 1) {
	  EraseFixedMutations(muttable[m].position, muttable[m].block, hh);
	  EraseFixedMutations(muttable[otherm_1].position, muttable[otherm_1].block, hh);
	  EraseFixedMutations(muttable[otherm_2].position, muttable[otherm_2].block, hh);

	  make_false_count += 1;

	  muttable[otherm_1].block = -9;
	  muttable[otherm_2].block = -9;
	} else if (muttable[m].frequency == 0 && muttable[otherm_1].frequency == 0 && muttable[otherm_2].frequency == 0) {
	  multihit[muttable[m].position] = false;
	  //	  multihit[muttable[otherm_1].position] = false;
	  //	  multihit[muttable[otherm_2].position] = false;

	  make_false_count += 1;
	  
	  muttable[otherm_1].block = -9;
	  muttable[otherm_2].block = -9;
	} else if (
		   (muttable[m].frequency >= 0 && muttable[m].frequency <= 1) ||
		   (muttable[otherm_1].frequency >= 0 && muttable[otherm_1].frequency <= 1) ||
		   (muttable[otherm_2].frequency >= 0 && muttable[otherm_2].frequency <= 1)
		    ) {
	  // SEGREGATING
	  temporalmuttable[mm].block = muttable[m].block;
	  temporalmuttable[mm].frequency = muttable[m].frequency;
	  temporalmuttable[mm].position = muttable[m].position;
	  mm++;
	  
	  temporalmuttable[mm].block = muttable[otherm_1].block;
	  temporalmuttable[mm].frequency = muttable[otherm_1].frequency;
	  temporalmuttable[mm].position = muttable[otherm_1].position;
	  mm++;

	  temporalmuttable[mm].block = muttable[otherm_2].block;
	  temporalmuttable[mm].frequency = muttable[otherm_2].frequency;
	  temporalmuttable[mm].position = muttable[otherm_2].position;
	  mm++;

	  muttable[otherm_1].block = -9;
	  muttable[otherm_2].block = -9;
	} else {
	  cout << muttable[m].frequency << " " << muttable[otherm_1].frequency << " " << muttable[otherm_2].frequency << '\n';
	  cout << "FSL_since_phaseIV() wrong here\n";
	  exit(0);
	}
      }
	// WHC: consider a senario with frequencies: 1, 0.5, 0 for m, otherm_1, otherm_2 !!!
	else {
	  cout << "something is wrong in FSL() here.\n";
	  exit(0);
	}

	/* WHC: this clause is a modificate of if (duFreq == true && duFreq_2 == false) */
	/* ================================================================ */
    }
  }
  MutCount = mm;
  for (m = 0; m < MutCount; m++) {
    muttable[m].block = temporalmuttable[m].block;
    muttable[m].frequency = temporalmuttable[m].frequency;
    muttable[m].position = temporalmuttable[m].position;
  }

  // WHC: just in case
  for (m = MutCount; m < MUTTABLESIZE; ++m) {
    muttable[m].block = -9;
    muttable[m].frequency = -9;
    muttable[m].position = -9;
  }


  /* TESTING */
  /*
  cout << "This time there are " << make_false_count << " multihit made to false.\n";

  cout << "The MutCount is " << MutCount << '\n';
  int num_of_true_in_multihit = 0, num_of_false_in_multihit = 0;
  for (int i = 0; i < BLOCKLENGTH; ++i) {
    if (multihit[i] == true) {
      num_of_true_in_multihit++;
    } else if (multihit[i] == false) {
      num_of_false_in_multihit++;
    } else {
      cout << "WRONG WRONG\n";
      exit(0);
    }
  }
  cout << "The number of true in multihit is " << num_of_true_in_multihit << " and false is " << num_of_false_in_multihit << endl;

  int num_of_true_in_multihit_single = 0, num_of_false_in_multihit_single = 0;
  for (int i = 0; i < BLOCKLENGTH; ++i) {
    if (multihit_single[i] == true) {
      num_of_true_in_multihit_single++;
    } else if (multihit_single[i] == false) {
      num_of_false_in_multihit_single++;
    } else {
      cout << "WRONG WRONG\n";
      exit(0);
    }
  }
  cout << "The number of true in multihit_single is " << num_of_true_in_multihit_single << " and false is " << num_of_false_in_multihit_single << endl;
  */
  
}


/////////////////////////////////////////////////////
////  COPYCHR, LOCATION & ERASE FIXED MUTATIONS  ////
/////////////////////////////////////////////////////

// COPY A CHROMOSOME (WHEN THERE IS NO RECOMBINATION)
void copychr(int prev, int ind0, int pres, int ind1) { // (origin,end)
  chrom *c1;
  int j, k;
  c1 = pointer[pres][ind1];
  c1->b = pointer[prev][ind0]->b;
  for (j = 0; j < pointer[prev][ind0]->b; j++) {
    c1->mpb[j] = pointer[prev][ind0]->mpb[j];
    for (k = 0; k < pointer[prev][ind0]->mpb[j]; k++) {
      c1->mutation[j][k] = pointer[prev][ind0]->mutation[j][k];
    }
  }

  // WHC: I think this is necessary, and I don't know if their original code is correct without this???
  //  if (c1->b == 2) {
  if (c1->b == 2) {
    c1->mpb[2] = 0;
    c1->mpb[3] = 0;
    c1->mpb[4] = 0;
    c1->mpb[5] = 0;
  } else if (c1->b == 4) {
    c1->mpb[4] = 0;
    c1->mpb[5] = 0;
  } else if (c1->b == 5) {
    cout << "hi, this should not be here! go for copychr_for_phaseVI()!\n";
    exit(0);
  } else if (c1->b == 6) {
    // do nothing
  } else {
    cout << "something wrong here.\n";
    exit(0);
  }
}

/* ================================================================ */
// WHC: a copychr() for phaseVI

void copychr_for_phaseVI(int prev, int ind0, int pres, int ind1) { // (origin,end)
  chrom *c1;
  int j, k;
  c1 = pointer[pres][ind1];
  c1->b = pointer[prev][ind0]->b;
  // WHC: similar as copychr(), this should work
  for (j = 0; j < 6; j++) {
    c1->mpb[j] = pointer[prev][ind0]->mpb[j];
    for (k = 0; k < pointer[prev][ind0]->mpb[j]; k++) {
      c1->mutation[j][k] = pointer[prev][ind0]->mutation[j][k];
    }
  }
}

/* ================================================================ */

// LOCATE A GIVEN MUTATION OR POSITION INSIDE A mutation VECTOR OF A GIVEN CHROMOSOME
int location(int position, int h, int ind, int j) {// Found the position in the c->mutation array where "position" should be
  chrom *c;
  int k = 0;
  c = pointer[h][ind];
  bool found = false;
  while (found == false && k <= c->mpb[j]) {
  // WHC: I think this is WRONG!!! ANd causes trouble; because it will return k = c->mpb[j], causes an undefined value!!
  // WHC: that's why void DivergenceForAll() could sometimes cause troubles (possibly)
  //  while (found == false && k < c->mpb[j]) {
  //  WHC: I don't know, should be careful before changing
    // WHC: I see; they double checked j < pionter[h][i]->mpb[blockA] in DivergenceForAll(), so that they know
    // WHC: thus, k == mpb[j] is a judgement; when k == mpb[j], it means "no found"
    if (position <= c->mutation[j][k] && k < c->mpb[j]) {
      found = true;
    }
    k++;
  }
  return (k - 1);
}

// ERASES MUTATIONS
void EraseFixedMutations(int fixed, int block, int prev) {
  int ii, k;

  if (block != 1) {
  multihit[fixed] = false;
  } else {
    multihit_single[fixed] = false;
  }
  
  for (ii = 0; ii < 2 * N; ii++) {
    if (pointer[prev][ii]->mpb[block] == 0) { continue; } // WHC: vital for this
    pointer[prev][ii]->mpb[block]--;
    for (k = location(fixed, prev, ii, block); k < pointer[prev][ii]->mpb[block]; k++) {
      pointer[prev][ii]->mutation[block][k] = pointer[prev][ii]->mutation[block][k + 1];
    }
  }
}


///////////////////////////////////////////////////////////////
////  SITE FREQUENCY SPECTRUM, BY BINS SFS, COLLAPSED SFS  ////
///////////////////////////////////////////////////////////////

// Calculates Site Frequency Spectrum with Printing option for mutationNew
// also prints SFS for blocks 0, 1 and 2
// also calculates pi, S, eta, H
float * SiteFrequencySpectrumPrint(int h, int block, int n, bool does_print) {
  int s = n;
  int p, i, j, k, index, m, number, duplicationFreq;
  //  int xi[2 * s], list [2 * N];
  int xi[2 * s], list [MutCount];
  // WHC: any good reason for list[] size less than 2*N??????
  float value, mufreq;
  static float results[] = {0, 0};

  results[0] = 0;
  results[1] = 0;

  if(block == 3){
    //duplicationFreq = Freq(prev, 2, n);
    duplicationFreq = DupliFreq(h, 3, n);
    // WHC: why divided by 2?? Because DupliFreq() returns the number of chroms, but s is the number of individuals!!
    s = (int) duplicationFreq/2;
  } else if (block == 5) {	// WHC: newly added
    duplicationFreq = DupliFreq(h, 5, n);
    s = (int) duplicationFreq / 2;
  }
  if (s == 0) { return results; }

  for (m = 0, number = 0; m < MutCount; m++) {
    if (block == muttable[m].block && muttable[m].frequency != 0) {
      if (muttable[m].frequency != 1) {
	// WHC: muttable[].frequency, is calculated in FSL() by dividing (2 * N), not (2 * s)
	// WHC: so, it is a population-level frequency, not a sample-level frequency as here
	// WHC: and this time, muFrequencyIntWholePopAndSample() will just do sample-level based on sample[]
	mufreq = (float) muFrequencyIntWholePopAndSample(h, muttable[m].position, muttable[m].block,n);
	mufreq = (mufreq)/(2*s);
	if(mufreq != 0){
	  list[number] = muttable[m].position;
	  number++; // number = segregatingSites in sample
	}
      }
    }
  }
  results[0] = number;
  // number=5;
  if(does_print == true){
    mutationsNewFile[block] << "\n\n//\nsegsites: "<< number <<"\npositions: ";
  }

  if (number == 0) { return results; }

  int prototype[number];
  int mutationCounts[number];
  for (index = 0; index < 2 * s; index++) {
    xi[index] = 0;
  }

  // BUILD THE PROTOTYPE: AN ORDERED SEQUENCE WITH ALL POSSIBLE MUTATION POSITIONS
  // WHC: sorting prototype[]
  for (m = 0, p = 0; m < number; m++) {
    if (p == 0) {
      prototype[p] = list[m];
      p++;
    }else {
      j = 0;
      while (list[m] > prototype[j] && j < p) { j++; }
      if (j < p) {
	for (k = p + 1; k > j; k--) {
	  prototype[k] = prototype[k - 1];
	}
      }
      prototype[j] = list[m];
      p++;
    }
  }

  // PRINTS POSITIONS OF SEGREGATING SITES WITH FIVE DECIMAL POSITIONS
  // WHC: so the output is ordered
  if(does_print == true){
    float rounded_mutation;
    for(int nn=0; nn < number ; nn++){
      rounded_mutation = (float) prototype[nn]/BLOCKLENGTH;
      rounded_mutation = round(rounded_mutation,7);
      mutationsNewFile[block] << rounded_mutation << " ";
    }
    mutationsNewFile[block] << "\n";
  }

  // THIS FUNCTION PRINTS THE CONTENT OF MUTATIONSNEWFILE
  if(does_print == true){
    int pos;
    for (i = 0; i < 2 * s; i++) {
      for (j = 0; j < number; j++) {
	pos = prototype[j];
	k = location(prototype[j], h, sample[i], block);
	if (pointer[h][sample[i]]->mutation[block][k] == pos && k < pointer[h][sample[i]]->mpb[block]) {
	  // WHC: this if () clause shows multiple times, it is to test whether the specific mutation (at pos position),
	  // WHC: shoulds in the specific individual (sample[i]) or not
	  // WHC: notice that the k < pointer[h][sample[i]]->mpb[block] is NECESSARY (hint: if there is no such mutation, return value..)
	  
	  mutationsNewFile[block] << "1";
	}else{
	  mutationsNewFile[block] << "0";
	}
      }
      mutationsNewFile[block] << "\n";
    }
  }

  //FIND THE FREQUENCY OF EACH MUTATION
  //MUTCOUNTS IS AN ARRAY THAT COUNTS THE NUMBER OF TIMES EACH MUTATION APPEARS
  for (j = 0; j < number; j++) {
    mutationCounts[j] = muFrequencyIntWholePopAndSample(h, prototype[j], block, n);
  }
  //FIND THE UNFOLDED SITE FREQUENCY SPECTRUM
  //xi[i] IS AN ARRAY THAT COUNTS THE NUMBER OF POLYMORPHIC SITES THAT HAVE i COPIES OF THE MUTATION
  // WHC: so the output file SFS_x_ID means -- 12, 10, 20,... there are 12 polymophic sites with 1 copy, 10 with 2 copies, 20 with 3...
  
  for (i = 1; i < 2 * s; i++) {
    for (j = 0; j < number; j++) {
      if (mutationCounts[j] == i) {
	xi[i]++;
      }
    }
  }

  if(does_print==true){
    for (i = 1; i < 2 * s; i++) {
      SFS[block] << xi[i] << " ";
    }
    SFS[block] << "\n";
  }

  // CALCULATE AVERAGE PAIRWISE DIFFERENCE
  // WHC: what does this mean??
  for (i = 1, value = 0; i < 2 * s; i++) {
    value += (float) (i * (2 * s - i) * xi[i]);
  }
  value /= (s * (2 * s - 1));
  results[1] = value;
  return results;

}

/* ================================================================ */
// WHC: new SiteFrequencySpectrumPrint() for phaseVI

float * SiteFrequencySpectrumPrint_for_phaseVI(int h, int block, int n, bool does_print) {
  int s = n;
  int p, i, j, k, index, m, number, duplicationFreq;
  //  int xi[2 * s], list [2 * N];
  // WHC: I think list[2*N] is incorrect
  int xi[2 * s], list[MutCount];
  float value, mufreq;
  static float results[] = {0, 0};

  results[0] = 0;
  results[1] = 0;

  if(block == 3){
    //duplicationFreq = DupliFreq(prev, 2, n);
    // WHC: duplicationFreq is the number of chroms that are carrying dup_1, 2 * N - duplicationFreq = #chroms without block 2 (dup_1)
    duplicationFreq = DupliFreq_for_phaseVI(h, 3, n);

    // WHC: why divided by 2?? Because DupliFreq() returns the number of chroms, but s is the number of individuals!!
    s = (int) duplicationFreq/2;
  } else if (block == 5) {	// WHC: newly added
    // WHC: as this is in phaseVI(), don't need this
    if (s * 2 != DupliFreq(h, 5, n)) { cout << "Hi, there is an error.\n"; exit(0); }
  }
  
  if (s == 0) {
    // WHC: as this is in pahseVI(), this should not happen (and also due to my limitation that chroms not carrying dup_1 ranges from 1 to (2 * N - 1)
    // WHC: wrong!! as we are sampling, not the whole population; so there could be possibilities when s == 0
    //    cout << "something wrong here here.\n";
    //    exit(0);
    return results;
  }

  for (m = 0, number = 0; m < MutCount; m++) {
    if (block == muttable[m].block && muttable[m].frequency != 0) {
      if (muttable[m].frequency != 1) {
	// WHC: muttable[].frequency, is calculated in FSL() by dividing (2 * N), not (2 * s)
	// WHC: so, it is a population-level frequency, not a sample-level frequency as here
	// WHC: and this time, muFrequencyIntWholePopAndSample() will just do sample-level based on sample[]
	mufreq = (float) muFrequencyIntWholePopAndSample(h, muttable[m].position, muttable[m].block,n);
	mufreq = (mufreq)/(2*s);
	if(mufreq != 0){
	  list[number] = muttable[m].position;
	  number++; // number = segregatingSites in sample
	}
      }
    }
  }
  results[0] = number;
  // number=5;
  if(does_print == true){
    mutationsNewFile[block] << "\n\n//\nsegsites: "<< number <<"\npositions: ";
  }

  if (number == 0) { return results; }

  int prototype[number];
  int mutationCounts[number];
  for (index = 0; index < 2 * s; index++) {
    xi[index] = 0;
  }

  // BUILD THE PROTOTYPE: AN ORDERED SEQUENCE WITH ALL POSSIBLE MUTATION POSITIONS
  // WHC: sorting prototype[]
  for (m = 0, p = 0; m < number; m++) {
    if (p == 0) {
      prototype[p] = list[m];
      p++;
    }else {
      j = 0;
      while (list[m] > prototype[j] && j < p) { j++; }
      if (j < p) {
	for (k = p + 1; k > j; k--) {
	  prototype[k] = prototype[k - 1];
	}
      }
      prototype[j] = list[m];
      p++;
    }
  }

  // PRINTS POSITIONS OF SEGREGATING SITES WITH FIVE DECIMAL POSITIONS
  // WHC: so the output is ordered
  if(does_print == true){
    float rounded_mutation;
    for(int nn=0; nn < number ; nn++){
      rounded_mutation = (float) prototype[nn]/BLOCKLENGTH;
      rounded_mutation = round(rounded_mutation,7);
      mutationsNewFile[block] << rounded_mutation << " ";
    }
    mutationsNewFile[block] << "\n";
  }

  // THIS FUNCTION PRINTS THE CONTENT OF MUTATIONSNEWFILE
  if(does_print == true){
    int pos;
    for (i = 0; i < 2 * s; i++) {
      for (j = 0; j < number; j++) {
	pos = prototype[j];
	k = location(prototype[j], h, sample[i], block);
	if (pointer[h][sample[i]]->mutation[block][k] == pos && k < pointer[h][sample[i]]->mpb[block]) {
	  // WHC: this if () clause shows multiple times, it is to test whether the specific mutation (at pos position),
	  // WHC: shoulds in the specific individual (sample[i]) or not
	  // WHC: notice that the k < pointer[h][sample[i]]->mpb[block] is NECESSARY (hint: if there is no such mutation, return value..)
	  
	  mutationsNewFile[block] << "1";
	}else{
	  mutationsNewFile[block] << "0";
	}
      }
      mutationsNewFile[block] << "\n";
    }
  }

  //FIND THE FREQUENCY OF EACH MUTATION
  //MUTCOUNTS IS AN ARRAY THAT COUNTS THE NUMBER OF TIMES EACH MUTATION APPEARS
  for (j = 0; j < number; j++) {
    mutationCounts[j] = muFrequencyIntWholePopAndSample(h, prototype[j], block, n);
  }
  //FIND THE UNFOLDED SITE FREQUENCY SPECTRUM
  //xi[i] IS AN ARRAY THAT COUNTS THE NUMBER OF POLYMORPHIC SITES THAT HAVE i COPIES OF THE MUTATION
  // WHC: so the output file SFS_x_ID means -- 12, 10, 20,... there are 12 polymophic sites with 1 copy, 10 with 2 copies, 20 with 3...
  
  for (i = 1; i < 2 * s; i++) {
    for (j = 0; j < number; j++) {
      if (mutationCounts[j] == i) {
	xi[i]++;
      }
    }
  }

  if(does_print==true){
    for (i = 1; i < 2 * s; i++) {
      SFS[block] << xi[i] << " ";
    }
    SFS[block] << "\n";
  }

  // CALCULATE AVERAGE PAIRWISE DIFFERENCE
  // WHC: what does this mean??
  for (i = 1, value = 0; i < 2 * s; i++) {
    value += (float) (i * (2 * s - i) * xi[i]);
  }
  value /= (s * (2 * s - 1));
  results[1] = value;
  return results;

}

/* ================================================================ */


////////////////////////////////////////////////////////////////
////  DIVERGENCE BETWEEN IN SAME AND DIFFERENT CHROMOSOMES  ////
////////////////////////////////////////////////////////////////

// CALCULATES DIVERGENCE BETWEEN TWO BLOCKS IN THE SAME CHROMOSOME
void DivergenceForAll(int blockA, int blockB, int h) {
  int i, k, j;
  int mutBlockA = 0;
  int mutBlockB = 0;
  int mutSharedAB = 0;
  int mutSharedABComp = 0;

  for (i = 0; i < 2 * N; i++) {
    for (k = 0; k < pointer[h][i]->mpb[blockA]; k++) {
      j = location(pointer[h][i]->mutation[blockA][k], h, i, blockB);
      if ((pointer[h][i]->mutation[blockA][k] == pointer[h][i]->mutation[blockB][j]) && (j < pointer[h][i]->mpb[blockB])) {
	mutSharedAB++;
      } else {
	mutBlockA++;
      }
    }
    for (k = 0; k < pointer[h][i]->mpb[blockB]; k++) {
      j = location(pointer[h][i]->mutation[blockB][k], h, i, blockA);
      if ((pointer[h][i]->mutation[blockB][k] == pointer[h][i]->mutation[blockA][j]) && (j < pointer[h][i]->mpb[blockA])) {
	mutSharedABComp++;
      } else {
	mutBlockB++;
      }
    }
  }
  if (mutSharedAB != mutSharedABComp) {
    cout << "ERROR en divergence for all\n" << mutSharedAB << " " << mutSharedABComp << "\n\n";
    //WHC: test if two pointer[][]->mutation[k or i] are the same
    for (int i = 0; i < 2 * N; ++i) {
      cout << "Block A (first 2):\n";
      cout << pointer[h][i]->mutation[blockA][0] << " " << pointer[h][i]->mutation[blockA][1] << '\n';
      //      for (int k = 0; k < pointer[h][i]->mpb[blockA]; ++k) {
      //	cout << pointer[h][i]->mutation[blockA][k] << " ";
      //      }
      cout << "Block B (first 2):\n";
      cout << pointer[h][i]->mutation[blockB][0] << " " << pointer[h][i]->mutation[blockB][1] << '\n';

      int k = pointer[h][i]->mpb[blockA] - 1;
      cout << "Block A (last 2):\n";
      cout << pointer[h][i]->mutation[blockA][k - 1] << " " << pointer[h][i]->mutation[blockA][k] << '\n';

      //      cout << '\n';
      cout << "Block B (last 2:)\n";

      k = pointer[h][i]->mpb[blockB] - 1;
      cout << pointer[h][i]->mutation[blockB][k - 1] << " " << pointer[h][i]->mutation[blockB][k] << '\n';

      cout << '\n';

    }
  }
}


////////////////////////////
////  SIMPLE FUNCTIONS  ////
////////////////////////////

// CALCULATES ABSOLUTE FREQUENCY OF THE DUPLICATION
// WHC: good for phaseIV(), after changing == to >=
int DupliFreq(int h, int block, int n) {
  int i = 0, quantity = 0;

  // WHC: Will this work for dup_2 after losing dup_1???

  if (loseFreq == false) {
    if (n == N){
      for (i = 0; i < 2 * n; i++) {
	//      if (pointer[h][i]->b == (block + 1)) {
	// WHC: maybe?
	if (pointer[h][i]->b >= (block + 1)) {
	  quantity++;
	}
      }
    } else {
      for (i = 0; i < 2 * n; i++) {
	//      if (pointer[h][sample[i]]->b == (block + 1)) {
	// WHC: same reason as previous
	if (pointer[h][sample[i]]->b >= (block + 1)) {
	  quantity++;
	}
      }
    }
  } else if (loseFreq == true) {
    // WHC: this only need to work for counting block = 4, since block = 2 is counted by DupliFreq_for_phaseVI()
    if (block != 5) { cout << "Does not do the job\n"; exit(0); }
    
    if (n == N){
      for (i = 0; i < 2 * n; i++) {
	//      if (pointer[h][i]->b == (block + 1)) {
	// WHC: maybe?
	if (pointer[h][i]->b >= 5) {
	  quantity++;
	}
      }
    } else {
      for (i = 0; i < 2 * n; i++) {
	//      if (pointer[h][sample[i]]->b == (block + 1)) {
	// WHC: same reason as previous
	if (pointer[h][sample[i]]->b >= 5) {
	  quantity++;
	}
      }
    }
  } else {
    cout << "no way to be here\n";
    exit(0);
  }
  
  return quantity;
}

/* ================================================================ */
// WHC: DupliFreq() for phaseVI

int DupliFreq_for_phaseVI(int h, int block, int n) {
  // WHC: can only work for counting dup_1 frequency in phaseVI
  int i = 0, quantity = 0;

  if (loseFreq != true || block != 3) { cout << "this boy does not do this job.\n"; exit(0); }
  
  if (n == N){
    for (i = 0; i < 2 * n; i++) {
      //      if (pointer[h][i]->b == (block + 1)) {
      // WHC: maybe?
      if (pointer[h][i]->b == 5) { // WHC: lost dup_1
	quantity++;
      }
    }
  }else{
    for (i = 0; i < 2 * n; i++) {
      //      if (pointer[h][sample[i]]->b == (block + 1)) {
      // WHC: same reason as previous
      if (pointer[h][sample[i]]->b == 5) { // WHC: lost dup_1
	quantity++;
      }
    }
  }
  return (2 * n - quantity);
}


/* ================================================================ */

// THIS FUNCTION SERVES TO SEARCH FOR THE ABSOLUTE FREQUENCY OF A GIVEN MUTATION IN THE WHOLE POPULATION OR A SAMPLE
int muFrequencyIntWholePopAndSample(int h, int pos, int block, int n) {
  int i, quantity, k;
  if (n == N) {
    for (i = 0, quantity = 0; i < 2 * n; i++) {
      k = location(pos, h, i, block);
      if (pointer[h][i]->mutation[block][k] == pos && k < pointer[h][i]->mpb[block]) {
	quantity++;
      }
    }
  } else {
    for (i = 0, quantity = 0; i < 2 * n; i++) {
      // WHC: imagine, how DISASTROUS it would be, if pointer[h]sample[i]] does NOT have block 4 (or, dup_2)!!!
      // WHC: the mpb[4] is undefined!!!
      // WHC: how did they overcome this for phaseII and block 2???
      // WHC: but they initiated this in int main()...?

      // WHC: I think start at phaseIV(), when a #blocks = 5 chrom in the table[4*N] has a "next-generation chrom" whose #blocks = 4
      // WHC: then, its mpb[4] is the same as previous, which is not initialized
      // WHC: which means, I FAILED to clear those values to 0!!!!!
      k = location(pos, h, sample[i], block);
      if (pointer[h][sample[i]]->mutation[block][k] == pos && k < pointer[h][sample[i]]->mpb[block]) {
	quantity++;
      }
    }
  }
  return quantity;
}

// THIS FUNCTION DOES THE CALLING FOR COLLAPSED DUPLICATIONS OF AN INDIVIDUAL FROM SAMPLE
int muFrequencyCollapsedCallingFromSample(int h, int pos, int n) {
  int i, j, quantity, collapsedQuantity, k, block;

  for (i = 0, collapsedQuantity=0; i < n; i++) {
    for(j = 0, quantity=0; j < 2; j++){
      block=0;
      k = location(pos, h, sample[(2*i)+j], block);
      if (pointer[h][sample[(2*i)+j]]->mutation[block][k] == pos && k < pointer[h][sample[(2*i)+j]]->mpb[block]) {
	quantity++;
      }
      block=2;
      if(pointer[h][sample[(2*i)+j]]->b > 2){
	k = location(pos, h, sample[(2*i)+j], block);
	if (pointer[h][sample[(2*i)+j]]->mutation[block][k] == pos && k < pointer[h][sample[(2*i)+j]]->mpb[block]) {
	  quantity++;
	}
      }
    }
    if(quantity == 4){ collapsedQuantity+=2;}
    else{
      if(quantity !=0){ collapsedQuantity++;}
    }
  }
  return collapsedQuantity;
}

//This function gives the quantity (integer) of a particular mutation found in a sample
//arguments are h (which chromosome array to search), pos (the integer indicating the mutation position),
//block (the block in which to search), n (sample size), init (the first individual chromosome to search
//in the sample), and jumps (the intervals with which to search within the sample)
int muFrequencySampleIntDiscont(int h, int pos, int block, int s, int init, int jumps) {
  int i, quantity, k;
	
  for (i = init, quantity = 0; i < 2 * s; i += jumps) {
    k = location(pos, h, sample[i], block);
    if (pointer[h][sample[i]]->mutation[block][k] == pos && k < pointer[h][sample[i]]->mpb[block]) {
      quantity++;
    }
  }
  return quantity;
}

// SEARCH A GIVEN MUTATION INSIDE MUTTABLE (RETURNS ITS POSITION)
int SearchMutation(int block, int pos, int mutcount) {
  int counter;
  bool found;

  counter = 0;
  found = false;
  while (found == false && counter < mutcount) {
    if (muttable[counter].position == pos && muttable[counter].block == block) {
      found = true;
    }
    counter++;
  }

  return (counter - 1);
}

// NUMERICAL FUNCTION: TRACTPQL (RETURN THE TRACT LENGTH)
int tractpql(float meanTL)
{
  // DISTRIBUTION TYPICAL PARAMETERS
  float q = (float) 1 / meanTL;
  float threshold = 0;
  float p;
  int x = 0;
  float power;

  p = rand() / ((float) RAND_MAX + 1);
  do {
    x++;
    power = ((float) x) - 0.5;
    threshold += q * pow(1 - q, power - 1);
  } while (p > threshold);

  return x;
}

// FUNCTION THAT CREATES A LIST OF n CHROMOSOMES TO SAMPLE. WORKS AWFUL WHEN n IS CLOSE TO 2*N
void SamplingIndividuals(int n) {
  int count, val;
  bool samplelist[2 * N];

  if (n != N) {
    for (count = 0; count < 2 * N; count++) {
      samplelist[count] = false;
    }
    count = 0;
    while (count < n) {
      val = (int) (rand() % (N));
      if (samplelist[val] == false) {
	samplelist[2 * val] = true;
	samplelist[2 * val + 1] = true;
	count++;
      }
    }
    for (val = 0, count = 0; val < 2 * N; val++) {
      if (samplelist[val] == true) {
	sample[count] = val;
	count++;
      }
    }
  } else {
    for (count = 0; count < 2 * N; count++) {
      sample[count] = count;
    }
  }
}

///////////////////////////////
////  FIXATION TRAJECTORY  ////
///////////////////////////////

int GenerateFixationTrajectory(int maxTime, int fixationTime) {
  int return_var;
  if (fixationTime != 0){
    // GENERATE A FIXATION TRAJECTORY WITH A CONSTANT INCREASE IN POPOULATION SIZE AND FIXED TO fixationTime
    // WHC: This trajectory is the absolute # of duplicated chromosomes in each generation
    int tt;
    float interval = (float) 2 * N / fixationTime;
    float fixationTrajectoryFloat[maxTime];

    fixationTrajectoryFloat[0] = 1;
    fixationTrajectory[0] = 1;
    for (tt = 1; tt < fixationTime; tt++) {
      fixationTrajectoryFloat[tt] = (float) fixationTrajectoryFloat[tt - 1] + interval;
      fixationTrajectory[tt] = (int) (fixationTrajectoryFloat[tt]);
    }
    for (tt = fixationTime; tt < STRUCTURED + 1; tt++) {
      fixationTrajectory[tt] = 2 * N;
    }
    return_var = fixationTime;
  }
  else if (fixationTime == 0){
    // THE FIXATION TRAJECTORY IS GENERATED THROUGH THE ALGORITHM PROPOSED BY KIMURA
    // IN PNAS 77 (1), 1980, USING A PSEUDO-SAMPLING METHOD.
    // THIS IS A MODIFIED VERSION THAT ALWAYS MAINTAINS A DIPLOID POPULATION WITH THE DUPLICATION
    int endTime=0, j, success, tt;
    float u, p, varP, psv;
    float fixationTrajectoryFloat[maxTime]; // Probability of fixation of the duplication in each generation (DupInd/2N)(relative freq)

	
    // Initialize duplication
    for (tt = 0; tt < 2; tt++) {
      fixationTrajectoryFloat[tt] = (float) 1 / (2 * N);
      fixationTrajectory[tt] = (int) (fixationTrajectoryFloat[tt]*2 * N);
      fixationTrajectory[tt] = 1;
    }
    success = 0;
    tt = 1;

    // Calculate trajectory
    do {
      u = rand() / ((float) RAND_MAX + 1);
      p = fixationTrajectoryFloat[tt];
      varP = p * (1 - p) / (2 * N);
      psv = sqrt(3 * varP)*(2 * u - 1);
      fixationTrajectoryFloat[tt + 1] = p + psv;
      tt++;
      fixationTrajectory[tt] = (int) (fixationTrajectoryFloat[tt]*2 * N);

      // MODIFICATION TO ALWAYS HAVE A DIPLOID POPULATION WITH THE DUPLICATION
      // If duplication is already fixed
      if (fixationTrajectoryFloat[tt] >= 1) {
	success = 1;
	endTime = tt;
	for (j = endTime; j < maxTime; j++) {
	  fixationTrajectory[j] = 2 * N;
	}
      }
      // If maxTime was not enough or duplication has been lost, repeat the generation of the trajectory (control to ensure fixation)
      if ((fixationTrajectory[tt] < 1) || (tt == maxTime)) {
	tt = 1;
	fixationTrajectoryFloat[tt] = (float) 1 / (2 * N);
      }
    } while (success == 0);
    return_var = endTime;
  }
  return return_var;
}


/* ================================================================================================ */
// WHC: generate a random number of individuals with 4 blocks, ranging from 1 to (2*N - 1)
// WHC: purely random for now, may introduce some more advanced function in the future

void Generate_phaseVI_trajectory(int mode) {

  /* ================================================================ */
  /*

    WHC: when I use  for (int i = 1; i < (PHASE_VI_LENGTH); ++i) {} clause, the last generation always yields 0 dup_1
    solved when I use for (int i = 1; i < (PHASE_VI_LENGTH) + 1; ++i) {} clause
    !! As in their phaseII(), GenerateFixationTrajectory(STRUCTURED + 1, ), there is a +1 there!

  
  phaseVI_trajectory[0] = 1;
  
  //  for (int i = 1; i < (PHASE_VI_LENGTH); ++i) {
  for (int i = 1; i < (PHASE_VI_LENGTH) + 1; ++i) {
    phaseVI_trajectory[i] = 1 + (rand() % (2 * N - 1)); // generating [1-(2N-1)]
  }

  for (int i = 0; i < 100; ++i) {
    cout << "phaseVI_trajectory is:\t";
    cout << phaseVI_trajectory[i] << '\n';
  }
  cout << endl;

  return;  
  WHC: will abandom this random method for now

  /* ================================================================ */  

  if (mode == 1) {
    phaseVI_trajectory[0] = 1;
    for (int i = 1; i < (PHASE_VI_LENGTH) + 1; ++i) {
      if (phaseVI_trajectory[i - 1] < (2 * N) * (1.0 / 4)) { // WHC: when #chrom < 1/4 total chroms
	phaseVI_trajectory[i] = phaseVI_trajectory[i - 1] + 1;
      } else if (phaseVI_trajectory[i - 1] > (2 * N) * (3.0 / 4)) { // WHC: when #chrom > 3/4 total chroms
	phaseVI_trajectory[i] = phaseVI_trajectory[i - 1] - 1;
      } else {			// WHC: between 1/4 and 3/4
	phaseVI_trajectory[i] = (int) ((2 * N) * (1.0 / 8) + (rand() % ((int) (2 * N * 6.0 / 8)))); // WHC: generating 2*N*[1/8-7/8]; VERY rough for now
      }
    }
  } else if (mode == 2) {
    phaseVI_trajectory[0] = 1;
    for (int i = 1; i < (PHASE_VI_LENGTH) + 1; ++i) {
      if (phaseVI_trajectory[i - 1] < N) {
	phaseVI_trajectory[i] = phaseVI_trajectory[i - 1] + 1;
      } else {
	phaseVI_trajectory[i] = (int) (2 * N * (1.0 / 2));
      }
    }

  }
}


/* ================================================================================================ */

void print_fertility(){
  ofstream out;
  out.open("fertility.out");
  for (std::vector<fertility_info>::iterator it = fertility_list.begin() ; it != fertility_list.end(); ++it){
    cout << (*it).x << " " << (*it).y << "\n";
  }
}

float round(float number_to_round, int decimal_places) { // pow() doesn't work with unsigned
  return int(number_to_round * pow(10.0, decimal_places) + .50001) /  pow(10.0, decimal_places);
}

int minim(int n1, int n2) {
  if (n1 > n2) { return n2;}
  else { return n1; }
}

bool sortx (fertility_info i,fertility_info j) { 
  int tt1 = (i).x;
  int tt2 = (j).x;
  return (tt1<tt2); 
}

bool sorty (fertility_info i,fertility_info j) { 
  int tt1 = (i).y;
  int tt2 = (j).y;
  return (tt1<tt2); 	
}


// WHC: randomly pick 0, 1 or 2, with 
int pick_a_pair(int mode) {
  if (mode == 1) {		// mode 1, giving dup_1 lower chance
    int p = rand() % 10;	// generate 0 - 9
    if (p == 0) {
      return 0;
    } else if (p == 9) {
      return 2;
    } else if (p > 0 && p < 9) {
      return 1;
    } else {
      cout << "wrong in pick_a_pair().\n";
      exit(0);
    }
  } else if (mode == 2) {	// mode 2, pick equally 0, 1, 2
    int p = rand() % 3;
    return p;
  }
  
}

// WHC: this number_of_nucleotide_diff() is for counting how many nucleotides are different between blockA in chrom i, and blockB in chrom j; for calculating pairwise-divergence between any two blocks

float number_of_nucleotide_diff(int blockA, int blockB, int chrom_1, int chrom_2, int h) {
  // only cares about a pair of blockA and blockB; iteration will be done in statistics() and statistics_for_phaseVI()
  // chrom_1 is the chrom that blockA belongs to; chrom_2 is the chrom that blockB belongs to
  
  int j;
  int mutBlockA = 0, mutBlockB = 0, mutSharedAB = 0, mutSharedABComp = 0;

  // WHC: roughly test if chrom_1 contains blockA and chrom_2 contains blockB
  if ( ! (pointer[h][chrom_1]->b >= blockA && pointer[h][chrom_2]->b >= blockB)) {
    cout << "probably something wrong here.\n";
    cout << pointer[h][chrom_1]->b << " " << pointer[h][chrom_2]->b << '\n';
    exit(0);
  }
  
  for (int k = 0; k < pointer[h][chrom_1]->mpb[blockA]; ++k) {
    // iterate through blockA in chrom_1
    j = location(pointer[h][chrom_1]->mutation[blockA][k], h, chrom_2, blockB);
    if ((pointer[h][chrom_1]->mutation[blockA][k] == pointer[h][chrom_2]->mutation[blockB][j]) && (j < pointer[h][chrom_2]->mpb[blockB])) {
      mutSharedAB++;
    } else {
      mutBlockA++;
    }
  }

  for (int k = 0; k < pointer[h][chrom_2]->mpb[blockB]; ++k) {
    // iterate through blockB in chrom_2
    j = location(pointer[h][chrom_2]->mutation[blockB][k], h, chrom_1, blockA);
    if ((pointer[h][chrom_2]->mutation[blockB][k] == pointer[h][chrom_1]->mutation[blockA][j]) && (j < pointer[h][chrom_1]->mpb[blockA])) {
	mutSharedABComp++;
      } else {
	mutBlockB++;
      }
  }

  if (mutSharedAB != mutSharedABComp) {
    cout << "ERROR in number_of_nucleotide_diff()\n" << mutSharedAB << " " << mutSharedABComp << '\n';
    exit(0);
  }

  float per_nucleotide_diff = ((float) ( mutBlockA + mutBlockB)) / ((float) BLOCKLENGTH);
  // WHC: IMPORTANT! which one should I return???? Should I divide by BLOCKLENGTH or not???
  
  //  cout << mutBlockA << " " << mutBlockB << '\n';
  // return per_nucleotide_diff;
  return (mutBlockA + mutBlockB);
}

float pairwise_divergence_between(int h, int n, int blockA, int blockB) {
  // n = sampleN[o]
    int num_of_chrom_carrying_blockA = 0, num_of_chrom_carrying_blockB = 0;
    float results_pairwise_divergence_between = 0.0;

    if (blockA == blockB) {
      cout << "well, this was not tested yet.\n";
      exit(0);
    }
    
    // double check the number of chroms carrying blockA or blockB
    int num_blockA = 0, num_blockB = 0;
    if (duFreq == true && duFreq_2 == false) {
      num_blockA = 2 * n, num_blockB = 0;
      for (int i = 0; i < 2 * n; ++i) {
	if (pointer[h][sample[i]]->b == 4) { ++num_blockB; }
      }
    }
    
    for (int r_1 = 0; r_1 < 2 * n; ++r_1) {
      if (duFreq == false) {
	cout << "this cannot be done without duplicated blocks!\n";
	exit(0);
	
      } else if (duFreq == true && duFreq_2 == false) {
	//	num_of_chrom_carrying_blockA = 0; should not add this here
	num_of_chrom_carrying_blockB = 0;
	if (blockA != 0 || blockB != 3) { cout << "this blockA should only be 0 here and blockB only 3.\n"; exit(0); }
	num_of_chrom_carrying_blockA++; // WHC: this blockA should only be 0

	for (int r_2 = 0; r_2 < 2 * n; ++r_2) {
	  if (pointer[h][sample[r_2]]->b == 4) { // has blockB = 3
	    num_of_chrom_carrying_blockB++;
	    results_pairwise_divergence_between += number_of_nucleotide_diff(blockA, blockB, sample[r_1], sample[r_2], h);
	    //	    cout << results_pairwise_divergence_between << " ";
	  }
	}
	
      } else if (duFreq_2 == true && loseFreq == false) {
	num_of_chrom_carrying_blockA = 0;
	num_of_chrom_carrying_blockB = 0;
	if (blockA == 0) {
	  num_of_chrom_carrying_blockA = 2 * n;
	} else if (blockA == 3) {
	  num_of_chrom_carrying_blockA = 2 * n;
	} else {
	  // to make sure that no need to skip blockA, as chrom_1 always contains blockA
	  cout << "blockA should always be the smaller block.\n";
	  exit(0);
	}
	
	if (blockB == 3) {
	  num_of_chrom_carrying_blockB = 2 * n;
	} else if (blockB == 5) {
	  //	  cout << "initial num_of_chrom_carrying_blockB is " << num_of_chrom_carrying_blockB << '\n';
	  for (int temp = 0; temp < 2 * n; ++temp) {
	    if (pointer[h][sample[temp]]->b == 6) { ++num_of_chrom_carrying_blockB; }
	  }
	}

	for (int r_2 = 0; r_2 < 2 * n; ++r_2) {
	  if (pointer[h][sample[r_2]]->b >= blockB + 1) { // has blockB

	    results_pairwise_divergence_between += number_of_nucleotide_diff(blockA, blockB, sample[r_1], sample[r_2], h);
	    //	    cout << results_pairwise_divergence_between << " ";
	  }
	}
	
      } else if (loseFreq == true) {
	num_of_chrom_carrying_blockA = 0;
	num_of_chrom_carrying_blockB = 0;

	//	cout << "should use statistics_for_phaseVI()!\n";
	//	exit(0);
	// WHC: now, block A should only be 0 or 5
	if (blockA == 0) {
	  num_of_chrom_carrying_blockA = 2 * n;
	} else if (blockA == 5) {
	  num_of_chrom_carrying_blockA = 2 * n;
	} else if (blockA == 3) {
	  cout << "this time blockA should only be 0 or 5.\n";
	  exit(0);
	} else {
	  cout << "should be something wrong here.\n";
	  exit(0);
	}

	if (blockB == 3) {
	  for (int temp = 0; temp < 2 * n; ++temp) {
	    if (pointer[h][sample[temp]]->b == 6) {
	      ++num_of_chrom_carrying_blockB;
	    }
	  }
	} else if (blockB == 5) {
	  num_of_chrom_carrying_blockB = 2 * n;
	} else if (blockB == 0) {
	  cout << "blockB cannot be 0";
	  exit(0);
	} else {
	  cout << "should be a mistake here.\n";
	  exit(0);
	}

	for (int r_2 = 0; r_2 < 2 * n; ++r_2) {
	  if (blockB == 3 && pointer[h][sample[r_2]]->b == 5) { continue; }
	  
	  results_pairwise_divergence_between += number_of_nucleotide_diff(blockA, blockB, sample[r_1], sample[r_2], h);
	    //	    cout << results_pairwise_divergence_between << " ";
	}
	
      } else {
	cout << "could not enter here.\n";
	exit(0);
      }
    }

    
    if (duFreq == true && duFreq_2 == false) {
      if (num_of_chrom_carrying_blockA != 2 * n) { cout << "WRONGWRONG\n"; exit(0); }
      if (num_of_chrom_carrying_blockB != num_blockB) { cout << "wrong number of chromB.\n"; exit(0); }
    }
    if (num_of_chrom_carrying_blockA != 2 * n) {
      cout << "The num_of_chrom_carrying_blockA should always be 2 * n instead of " << num_of_chrom_carrying_blockA << '\n';
      exit(0);
    }
    
    if (num_of_chrom_carrying_blockA > 2 * n || num_of_chrom_carrying_blockB > 2 * n) {
      cout << num_of_chrom_carrying_blockA << " " << num_of_chrom_carrying_blockB << '\n';
      cout << "something wrong in counting.\n";
      exit(0);
    }

    if (num_of_chrom_carrying_blockA == 0 || num_of_chrom_carrying_blockB == 0) {
      // WHC: avoid any 0 numbers
      return 0.0;
    }
    
    results_pairwise_divergence_between = results_pairwise_divergence_between / (num_of_chrom_carrying_blockA * num_of_chrom_carrying_blockB);
    
    return results_pairwise_divergence_between;

}