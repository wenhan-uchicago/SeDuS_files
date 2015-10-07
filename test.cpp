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

int main() {



  // WHC: donorRatio initialization


    
const float ori_to_dup_1 = 0.5, ori_to_dup_2 = 0.5, dup_1_to_dup_2 = 0.5;
// Percentage of gene conversion events that occur from the original to the duplicated block
float donorRatio[5][5] = {
  {-1, -1, ori_to_dup_1, -1, ori_to_dup_2},
  {-1, -1, -1, -1, -1},
  {1 - ori_to_dup_1, -1, -1, -1, dup_1_to_dup_2},
  {-1, -1, -1, -1, -1},
  {1 - ori_to_dup_2, -1, 1 - dup_1_to_dup_2, -1, -1}
};
  float (*p)[5] = donorRatio;  
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      cout << donorRatio[i][j] << " ";
    }
    cout << '\n';
  }
  cout << p[4][2] << endl;

  return 0;
}
