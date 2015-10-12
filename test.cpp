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

int phaseVI_trajectory[1000];

void Generate_phaseVI_trajectory() {
  for (int i = 0; i < 1000; ++i) {
    phaseVI_trajectory[i] = 1 + (rand() % (2 * N - 1)); // generating [1-(2N-1)]
  }

  //  for (int i = 0; i < 1000; ++i) {
  //    cout << phaseVI_trajectory[i] << " ";
  //  }
  //  cout << endl;

}
