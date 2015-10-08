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

struct chrom {
  int b;
};

const float ori_to_dup_1 = 0.5, ori_to_dup_2 = 0.5, dup_1_to_dup_2 = 0.5;

float p_donorRatio[5][5] = {
  {-1, -1, ori_to_dup_1, -1, ori_to_dup_2},
  {-1, -1, -1, -1, -1},
  {1 - ori_to_dup_1, -1, -1, -1, dup_1_to_dup_2},
  {-1, -1, -1, -1, -1},
  {1 - ori_to_dup_2, -1, 1 - dup_1_to_dup_2, -1, -1}
};


int main() {

  cout << "donor\treceptor\tchrDonor\tchrReceptor\n";

  time_t seconds;	// WHC: seconds between now and 1970/1/1
  time(&seconds);	// WHC: get the abovementioned seconds
  srand((unsigned int) seconds);

  
  int i = 0, partner = 1;
  int block_1, block_2, donor, receptor, chrDonor, chrReceptor;


  float donorRatio;
  int ori_index = 0, dup_1 = 2, dup_2 = 4;
  chrom chr_1, chr_2;
  chrom *chr1 = &chr_1, *chr2 = &chr_2;
  chr1->b = 5;
  chr2->b = 3;
  float p;

  for (int m = 0; m < 10; ++m) {
      
  if (chr1->b == 3 && chr2->b <= 3) {		// WHC: pick 2 pairs
    block_1 = ori_index;
    block_2 = dup_1;
    donorRatio = p_donorRatio[0][2];
  } else if (chr1->b == 5 || chr2->b == 5) {	// WHC: could pick ori, dup_1, dup_2 blocks
    // WHC: randomly generate 0, 1, 2
    // WHC: which correspond to ori&dup_1, ori&dup_2 and dup_1&dup_2 pairs
    int which_pair = rand() % 3;
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
      
  if (chr2->b == 5 && chr1->b == 5) {
    p = rand() / ((float) RAND_MAX + 1); // WHC: randomly choose from which chromosome to which
    if (p < 0.5) {
      chrDonor = i;
      chrReceptor = partner;
    } else {
      chrDonor = partner;
      chrReceptor = i;
    }
  } else if (chr2->b == 5 && chr1->b == 3) {
    if (block_2 == dup_2) {	// only chr2 can be donor
      chrDonor = partner;
      chrReceptor = i;
    } else {		// randomly chosen
      p = rand() / ((float) RAND_MAX + 1); // WHC: randomly choose from which chromosome to which
      if (p < 0.5) {
	chrDonor = i;
	chrReceptor = partner;
      } else {
	chrDonor = partner;
	chrReceptor = i;
      }
    }
  } else if (chr2->b == 3) { // WHC: chr2 only contains ori and dup_1 (phase which loses dup_1 will not be considered here
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
    if (donor == dup_1 || donor == dup_2) { // this chrom does not have duplications
      chrDonor = i;
      chrReceptor = partner;
    } else {
      chrDonor = partner;
      chrReceptor = i;
    }
  }

  cout << donor << '\t' << receptor << '\t' << chrDonor << '\t' << chrReceptor << '\n';
}

  return 0;
}
