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
float * see_return_local_array();

int main() {

  float * p = see_return_local_array();
  return 0;
}

float * see_return_local_array() {
  static float result[] = {0, 0};
  return result;
}
