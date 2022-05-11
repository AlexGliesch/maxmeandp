
/*****************************************************************************\
 * This program aims to generate the benchmark instances used in our paper   *
 * for Max-mean Dispersion problem.                                          *
 * Modified by AlexGliesch on 2022-05-11 for research purposes.
\*****************************************************************************/

/*****************************************************************************/
/**********          0. Header files and variables      **********************/
/*****************************************************************************/
// #include <conio.h>
#include <ctime>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
using namespace std;
FILE* fp;
int N;      // node number in graph
double** D; // distance matrix

/*****************************************************************************/
/*****************          1. Outputing  results      ***********************/
/*****************************************************************************/
void Outputing(int tt) {
  int i, j;
  FILE* fp;
  char buff[80];

  sprintf(buff, "%s%d_%d.txt", "MDPI", tt, N);
  fp = fopen(buff, "a+");
  // fprintf(fp,"%d\n", N);
  for (i = 0; i < N; i++)
    for (j = i + 1; j < N; j++)
      fprintf(fp, "%d   %d   %5.2lf \n", i + 1, j + 1, D[i][j]);
  fclose(fp);
}

/*****************************************************************************/
/*****************         2. Main Function           ************************/
/*****************************************************************************/
int main(int argc, char** argv) {
  int i, j, i1, j1, seed;
  int tt;
  N = atoi(argv[1]); // 5000;  //N = 3000,5000.

  seed = (N == 3000 ? 1000 : 10000); // seed = 1000 for n=3000, seed = 10000 for n=5000
  // seed = 5000; // seed = 3000 for N =3000, seed = 5000 for N = 5000
  srand(seed);

  int sign;
  // N = 5000; // N = 3000, 5000.
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
        D[i][j] = sign * ((rand() % 30000) / 30000.0) * 10.0;
        D[j][i] = D[i][j];
      }
    Outputing(tt);
  }

  return 0;
}
