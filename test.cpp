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
int N = 100;

int main() {

    phaseVI_trajectory[0] = 1;
    for (int i = 1; i < 1000; ++i) {
      if (phaseVI_trajectory[i - 1] < (2 * N) * (2.0 / 4)) {
	phaseVI_trajectory[i] = phaseVI_trajectory[i - 1] + 1;
      } else {
	phaseVI_trajectory[i] = (int) (2 * N * (1.0 / 2));
      }
    }


  for (int i = 0; i < 1000; ++i) {
    cout << phaseVI_trajectory[i] << " ";
  }
  cout << endl;
  return 0;
}

