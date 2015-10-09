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
  chrom chr_1, chr_2;
  chrom *chr1 = &chr_2, *chr2 = &chr_2;
  chr1->b = 5;
  chr2->b = 3;
  float p;

  return 0;
}
