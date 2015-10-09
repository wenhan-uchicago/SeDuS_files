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

int BLOCKLENGTH = 10000; // Block length
#define maxNumOfHS 3  // Maximum number of crossover hotspots
// WHC: hard cord to initiate these values
int crossoverBegin[maxNumOfHS] = {0, 3 * BLOCKLENGTH, 4 * BLOCKLENGTH}; // Start point of crossover regions
int crossoverEnd[maxNumOfHS] = {3 * BLOCKLENGTH, 4 * BLOCKLENGTH, 5 * BLOCKLENGTH}; // End point of crossover regions
// WHC: artificially selected values, to represent 2MB region as crossover hotspot
float crossoverRatio[maxNumOfHS] = {0.15, 0.8, 0.05}; // Relative weights of crossover regions
int numHS = 3; // WHC: Number of hotspots, single_2 is a scaled hopspot representing a 2MB region

int main() {
  float total = 0;
  for (int HS = 0; HS < numHS; HS++) {
    cout << crossoverRatio[HS] << '\n';
    total += crossoverRatio[HS];
  }
  cout << total << endl;
  return 0;
}
