
/*****************************************************************************\
 * This program aims to generate the benchmark instances used in our paper   *
 * for the weighted Max-Mean dispersion problem.                             *
 * Modified by AlexGliesch on 2022-05-11 for research purposes.
\*****************************************************************************/

/*****************************************************************************/
/**********          0. Header files and variables      **********************/
/*****************************************************************************/
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
// #include <conio.h>
#include <ctype.h>
using namespace std;
FILE* fp;
int N;      // node number in graph
double** D; // distance matrix
double* w;  // node weight
const double WM = 5.0;
constexpr bool print_n = false;
constexpr bool print_weight = false;
/*****************************************************************************/
/*****************          1. Outputing  results      ***********************/
/*****************************************************************************/
void Outputing(int tt) {
  int i, j;
  FILE* fp;
  char buff[80];

  sprintf(buff, "%s_%d_%d.txt", "II", N, tt);
  fp = fopen(buff, "a+");
  if (print_n) fprintf(fp, "%d\n", N);
  if (print_weight)
    for (i = 0; i < N; i++)
      fprintf(fp, "%d   %5.2lf\n", i + 1, w[i]);
  for (i = 0; i < N; i++)
    for (j = i + 1; j < N; j++)
      fprintf(fp, "%d   %d   %5.2lf\n", i + 1, j + 1, D[i][j]);
  fclose(fp);
}

/*****************************************************************************/
/*****************         2. Main Function           ************************/
/*****************************************************************************/
int main(int argc, char** argv) {
  int i, j, i1, j1, seed;
  int tt;

  N = atoi(argv[1]);
  seed = N;
  //  seed = 5000 ; // seed = 1000 for n=1000, seed = 3000 for n=3000, seed = 5000 for
  //  n=5000
  srand(seed);

  int sign;
  //  N = 5000;     // N = 1000, 3000, 5000.
  w = new double[N];
  D = new double*[N];
  for (i = 0; i < N; i++)
    D[i] = new double[N];
  for (tt = 1; tt <= 10; tt++) {
    for (i = 0; i < N; i++)
      for (j = i + 1; j < N; j++) {
        if (rand() % 2)
          sign = 1;
        else
          sign = -1;
        D[i][j] = sign * (((rand() % 30000) / 30000.0) * 5.0 + 5.0);
        D[j][i] = D[i][j];
      }
    for (i = 0; i < N; i++)
      w[i] = 1 + WM * (rand() % 30000) / 30000.0;

    Outputing(tt);
  }

  return 0;
}
