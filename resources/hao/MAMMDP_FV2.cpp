
/*****************************************************************************\
 *                            Algorithm:  MA-MMDP                            *
 * This program demonstrates the use of memetic algorithm for solving        *
 * the Max-mean Diversity Problem (MMDP). For more information about         *
 * this algorithm (named as MA-MMDP), email to: laixiangjing@gmail.com       *
 *
 *
 * Modified by AlexGliesch on 2022-04-29 for research purposes.
\*****************************************************************************/

/*****************************************************************************/
/**********  0. Header files and global varialbes       **********************/
/*****************************************************************************/
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
const int Npop = 10; // size of population
char* File_Name;
FILE* fp;
char* outfilename;

double starting_time, total_time_suc, time_limit, time_one_run;
int Nhit;
int N, M; // node number in graph
double f, f_best, BestResult, AvgResult, N_times;
int* S; // common solution array
int* S_best;
int* SA;
int* SB;

double* P;               // the potential array about the set M
double** D;              // distance matrix
int* TabuTenure;         // tabu tenure
const int alpha = 50000; // the depth of tabu search

typedef struct Solution {
  int* S;
  double value;
} Solution;

typedef struct PairSet {
  int i;
  int j;
  int number;
} Pair;

Solution* POP;
Solution S_off;
Solution SS_BEST;
Solution G_BEST;

int** Pair_Set;
Pair* pair_s;

/*****************************************************************************/
/*****************          1. Initializing           ************************/
/*****************************************************************************/
// 1.1 Inputing
void Initializing() {
  int i, x1, x2;
  double d;

  ifstream FIC;
  FIC.open(File_Name);
  if (FIC.fail()) {
    cout << "### Erreur open, File_Name " << File_Name << endl;
    exit(0);
  }
  N = 0;

  while (!FIC.eof()) {
    FIC >> x1 >> x2 >> d;
    N = std::max(x1, x2);
  }
  FIC.clear();
  FIC.seekg(0);

  // FIC >> N;

  D = new double*[N];
  for (i = 0; i < N; i++)
    D[i] = new double[N];

  S = new int[N];
  S_best = new int[N];
  SA = new int[N];
  SB = new int[N];
  P = new double[N];
  TabuTenure = new int[N];

  while (!FIC.eof()) {
    FIC >> x1 >> x2 >> d;

    if (x1 < 0 || x2 < 0) {
      cout << "### Error of node : x1=" << x1 << ", x2=" << x2 << endl;
      exit(0);
    }
    D[x1 - 1][x2 - 1] = D[x2 - 1][x1 - 1] = d;
  }
  FIC.close();
}

void AssignMemery() {
  int i, j;
  POP = new Solution[Npop];
  for (i = 0; i < Npop; i++)
    POP[i].S = new int[N];
  S_off.S = new int[N];
  SS_BEST.S = new int[N];
  G_BEST.S = new int[N];
  Pair_Set = new int*[Npop];
  for (i = 0; i < Npop; i++)
    Pair_Set[i] = new int[Npop];
  for (i = 0; i < Npop; i++)
    for (j = 0; j < Npop; j++)
      Pair_Set[i][j] = 0;
  pair_s = new Pair[Npop * Npop / 2];
}

/*****************************************************************************/
/*****************  2. One-Flip Neighborhood  Tabu Search       **************/
/*****************************************************************************/
// 2.1 Clear the delta matrix
void Clear_Potential() {
  int x;
  f = 0.0;
  f_best = -99999;
  for (x = 0; x < N; x++)
    P[x] = 0.0;
  for (x = 0; x < N; x++)
    TabuTenure[x] = 0;
  return;
}
// 2.2 Build delta matrix
void Build_Potential() {
  int i, j;
  Clear_Potential();
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      if ((i != j) && (1 == S[j])) {
        P[i] += D[i][j];
      }
  M = 0;
  for (i = 0; i < N; i++)
    if (S[i] == 1) M++;
  f = 0;
  for (i = 0; i < N; i++)
    for (j = i + 1; j < N; j++)
      if ((S[i] == 1) && (S[j] == 1)) f += D[i][j];
  f = f / M;
  // printf("\n N =%d f=%lf M =%d \n",N,f,M);
  return;
}

// 2.3 Update potential vector after one-flip
void One_Flip_Update_Potential(int i, int old_value) {
  int j;
  for (j = 0; j < N; j++)
    if (j != i) {
      if (old_value == 0)
        P[j] += D[i][j];
      else
        P[j] -= D[i][j];
    }
  if (old_value == 0)
    M++;
  else
    M--;
  return;
}
// 2.4 Taboo Search with the Flip neighborhood N1
double One_Flip_Tabu_Search(int Sol[], double* value, int alpha) {
  int i, iter, num;
  int non_improve = 0;         // the stop condition of TS
  int num_tabu_best, num_best; // the number of tabu neighbors and non-tabu neighbors
  int best_x[50], x;
  int tabu_best_x[50];
  int select, old_value;
  double tabu_best_delta, best_delta, delt;

  int t, p;
  const int p_max = 15;
  int B[p_max] = {1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1};
  int A[p_max];
  int T[p_max];
  int T_max = 120;

  // double total_time, starting_time;
  // starting_time = clock();

  for (i = 0; i < p_max; i++) {
    A[i] = 5 * T_max * B[i] / 8;
    T[i] = T_max * B[i] / 8;
  }

  for (i = 0; i < N; i++)
    S[i] = Sol[i];
  Build_Potential();
  f_best = f;
  for (i = 0; i < N; i++)
    S_best[i] = S[i];

  // printf("\n\n");
  // cout << endl << "One_Move :       iter         f           f_best         time" <<
  // endl; cout << "---------------------------------------------------" << endl;

  p = 0;
  t = 0;
  iter = 0;
  while (non_improve < alpha) {
    num = 0;
    tabu_best_delta = -9999999.0;
    best_delta = -9999999.0;
    num_tabu_best = 0;
    num_best = 0;

    for (x = 0; x < N; x++) {
      num++;
      if (1 == S[x])
        delt = (-1.0 * P[x]) / (M - 1) + f / (M - 1); // To compute the vaule of delt.
      else
        delt = P[x] / (M + 1) - f / (M + 1);

      if (TabuTenure[x] <= iter) // if this is not tabued
      {
        if (delt > best_delta) {
          best_x[0] = x;
          best_delta = delt;
          num_best = 1;
        } else if (delt == best_delta && num_best < 50) {
          best_x[num_best] = x;
          num_best++;
        }
      } else if (TabuTenure[x] > iter) // if it is tabu
      {
        if (delt > tabu_best_delta) {
          tabu_best_x[0] = x;
          tabu_best_delta = delt;
          num_tabu_best = 1;
        } else if (delt == tabu_best_delta && num_tabu_best < 50) {
          tabu_best_x[num_tabu_best] = x;
          num_tabu_best++;
        }
      }
    }

    if ((num_tabu_best > 0 && tabu_best_delta > best_delta &&
         (f + tabu_best_delta > f_best)) ||
        num_best == 0) // aspiration criterion
    {
      f += tabu_best_delta;
      select = rand() % num_tabu_best;

      old_value = S[tabu_best_x[select]];
      One_Flip_Update_Potential(tabu_best_x[select], old_value);
      S[tabu_best_x[select]] = 1 - S[tabu_best_x[select]];

      TabuTenure[tabu_best_x[select]] = T[p] + rand() % 3;
      TabuTenure[tabu_best_x[select]] += iter;

      t++;
      if (t > A[p]) {
        p = (p + 1) % p_max;
        t = 0;
      }
    } else // non tabu
    {
      f += best_delta;
      select = rand() % num_best;

      old_value = S[best_x[select]];
      One_Flip_Update_Potential(best_x[select], old_value);
      S[best_x[select]] = 1 - S[best_x[select]];

      TabuTenure[best_x[select]] = T[p] + rand() % 3;
      TabuTenure[best_x[select]] += iter;

      t++;
      if (t > A[p]) {
        p = (p + 1) % p_max;
        t = 0;
      }
    }
    iter++;
    // total_time = (double) (1.0*(clock()-starting_time)/CLOCKS_PER_SEC);
    if (f > f_best + 10e-7) {
      // if( f > f_best + 10e-7 )
      //   {
      f_best = f;
      for (i = 0; i < N; i++)
        S_best[i] = S[i];
      // printf("\n One_Flip :  %8d       %lf       %lf     %lf", iter, f,
      // f_best,total_time);
      non_improve = 0;
      //   }
      // else if ( f == f_best )  non_improve ++ ;
    } else
      non_improve++;
  }
  for (i = 0; i < N; i++)
    Sol[i] = S_best[i];
  *value = f_best;
  return f_best;
}
/*****************************************************************************/
/*****************      3. Population Initialization  ************************/
/*****************************************************************************/
void InitiaSol(int* Sol) {
  int i;
  for (i = 0; i < N; i++)
    Sol[i] = rand() % 2;
} // Initial solution

void InitiaPoP(Solution* PP) {
  int i, j, s;
  for (i = 0; i < Npop; i++) {
    InitiaSol(PP[i].S);
    One_Flip_Tabu_Search(PP[i].S, &PP[i].value, alpha);
    if (PP[i].value > SS_BEST.value) {
      for (j = 0; j < N; j++)
        SS_BEST.S[j] = PP[i].S[j];
      SS_BEST.value = PP[i].value;
    }

    if (BestResult < 0.0) {
      BestResult = PP[i].value;
    } else if (PP[i].value > BestResult + 10e-7) {
      BestResult = PP[i].value;
      N_times = 0;
      total_time_suc = (double)((clock() - starting_time) / CLOCKS_PER_SEC);
    } else if (fabs(PP[i].value - BestResult) < 10e-7) {
      if (N_times == 0)
        total_time_suc += (double)((clock() - starting_time) / CLOCKS_PER_SEC);
      N_times++;
    }

    time_one_run = (double)((clock() - starting_time) / CLOCKS_PER_SEC);
    if (time_one_run > time_limit) return;
  }

  for (i = 0; i < Npop; i++)
    for (j = i + 1; j < Npop; j++)
      Pair_Set[i][j] = 1;
  s = 0;
  for (i = 0; i < Npop; i++)
    for (j = i + 1; j < Npop; j++)
      if (Pair_Set[i][j] == 1) {
        pair_s[s].i = i;
        pair_s[s].j = j;
        s++;
      }
  pair_s[0].number = s;
  // printf("intital number of Pairset : %d \n",pair_s[0].number);
}
void Rebuild(Solution* PP, Solution S_B) {
  int i, j, s;
  double worst_value = 999999999;
  int k;
  for (i = 0; i < Npop; i++) {
    InitiaSol(PP[i].S);
    One_Flip_Tabu_Search(PP[i].S, &PP[i].value, alpha);

    if (PP[i].value > BestResult + 10e-7) {
      BestResult = PP[i].value;
      N_times = 0;
      total_time_suc = (double)((clock() - starting_time) / CLOCKS_PER_SEC);
    } else if (fabs(PP[i].value - BestResult) < 10e-7) {
      if (N_times == 0)
        total_time_suc += (double)((clock() - starting_time) / CLOCKS_PER_SEC);
      N_times++;
    }

    if (PP[i].value < worst_value) {
      worst_value = PP[i].value;
      k = i;
    }
    time_one_run = (double)((clock() - starting_time) / CLOCKS_PER_SEC);
    if (time_one_run > time_limit) return;
  }
  for (j = 0; j < N; j++)
    PP[k].S[j] = S_B.S[j];
  PP[k].value = S_B.value;

  for (i = 0; i < Npop; i++)
    for (j = i + 1; j < Npop; j++)
      Pair_Set[i][j] = 1;
  s = 0;
  for (i = 0; i < Npop; i++)
    for (j = i + 1; j < Npop; j++)
      if (Pair_Set[i][j] == 1) {
        pair_s[s].i = i;
        pair_s[s].j = j;
        s++;
      }
  pair_s[0].number = s;
  // printf("rebuild intital number of Pairset : %d \n",pair_s[0].number);
} // rebuid the population
/*****************************************************************************/
/*****************          4. Outputing  results      ***********************/
/*****************************************************************************/

void OutSol(Solution& S, char* filename) {
  int i;
  int r;
  FILE* fp;
  char buff[80];
  sprintf(buff, "%s.sol", filename);
  fp = fopen(buff, "a+");
  fprintf(fp, "N= %d  f= %lf\n", N, S.value);
  for (i = 0; i < N; i++)
    fprintf(fp, "%d   %d\n", i + 1, S.S[i]);
  fclose(fp);
} // output a solution
void Outresulting(char* filename, char* outfile, double time) {
  int i, j;
  FILE* fp;
  char buff[80];
  sprintf(buff, "%s", outfile);
  fp = fopen(buff, "a+");
  fprintf(fp, "%s  %d   %lf   %lf  %d   %lf \n", filename, N, BestResult, AvgResult, Nhit,
          time);
  fclose(fp);
} // output the statistical result

/*****************************************************************************/
/*****************         5. Crossover Operator       ***********************/
/*****************************************************************************/
// 1) The greedy crossover operator
void Crossover(Solution S1, Solution S2, Solution S3) {
  int i, j;
  int i_index;
  int ma = 0;
  int mb = 0;
  int m = 0, mc;
  double delt, max_delt;

  for (i = 0; i < N; i++) {
    SA[i] = S1.S[i];
    SB[i] = S2.S[i];
    S3.S[i] = 0;
  }
  for (i = 0; i < N; i++) {
    if (SA[i] == 1) ma++;
    if (SB[i] == 1) mb++;
  }
  mc = (ma + mb) / 2;

  // printf("\n f1 =%lf  f2=%lf mc =%d  ma =%d  mb =%d \n",S1.value,S2.value, mc, ma, mb);
  for (i = 0; i < N; i++)
    if ((SA[i] == 1) && (SB[i] == 1)) {
      m++;
      S3.S[i] = 1;
      SA[i] = 0;
      ma--;
      SB[i] = 0;
      mb--;
    }

  // printf("\n common elements m =%d  ma =%d  mb =%d \n", m, ma, mb);
  for (i = 0; i < N; i++)
    P[i] = 0.0;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      if ((i != j) && (1 == S3.S[j])) {
        P[i] += D[i][j];
      }

  f = 0;
  for (i = 0; i < N; i++)
    for (j = i + 1; j < N; j++)
      if ((S3.S[i] == 1) && (S3.S[j] == 1)) f += D[i][j];
  f = f / m;

  while (m < mc) {
    if (ma > 0) {
      max_delt = -99999;
      for (i = 0; i < N; i++) {
        if (SA[i] == 1) {
          delt = P[i] / (m + 1) - f / (m + 1);
          if (delt > max_delt) {
            max_delt = delt;
            i_index = i;
          }
        }
      }
      S3.S[i_index] = 1;
      m++;
      SA[i_index] = 0;
      ma--;
      f += max_delt;
      for (j = 0; j < N; j++)
        if (j != i_index) P[j] += D[i_index][j];
    }

    if (mb > 0) {
      max_delt = -99999;
      for (i = 0; i < N; i++) {
        if (SB[i] == 1) {
          delt = P[i] / (m + 1) - f / (m + 1);
          if (delt > max_delt) {
            max_delt = delt;
            i_index = i;
          }
        }
      }
      S3.S[i_index] = 1;
      m++;
      SB[i_index] = 0;
      mb--;
      f += max_delt;
      for (j = 0; j < N; j++)
        if (j != i_index) P[j] += D[i_index][j];
    }
  }
  // printf("\n  f = %lf   m = %d\n",f,m);
}
// 2) The random crossover operator
void RandomCrossover(Solution S1, Solution S2, Solution S3) {
  int i, j, r;
  int i_index;
  int ma = 0, mb = 0, mc = 0;

  for (i = 0; i < N; i++) {
    SA[i] = S1.S[i];
    SB[i] = S2.S[i];
  }

  for (i = 0; i < N; i++) {
    r = rand() % 2;
    if (r == 1)
      S3.S[i] = S1.S[i];
    else
      S3.S[i] = S2.S[i];
  }
}
/***************************************************************************/
/****************************** 6. Updating PairSet ************************/
/***************************************************************************/
void updating_PairSet(int x, int y) {
  int i, j;
  int s = 0;
  Pair_Set[x][y] = 0;
  for (i = 0; i < Npop; i++)
    for (j = i + 1; j < Npop; j++)
      if (Pair_Set[i][j] == 1) {
        pair_s[s].i = i;
        pair_s[s].j = j;
        s++;
      }

  pair_s[0].number = s;
}

void updating_PairSet_pop(int position) {
  int i, j, s = 0;
  for (i = position + 1; i < Npop; i++)
    Pair_Set[position][i] = 1;
  for (j = 0; j < position; j++)
    Pair_Set[j][position] = 1;

  for (i = 0; i < Npop; i++)
    for (j = i + 1; j < Npop; j++)
      if (Pair_Set[i][j] == 1) {
        pair_s[s].i = i;
        pair_s[s].j = j;
        s++;
      }
  pair_s[0].number = s;
}
/*****************************************************************************/
/*****************         7. Pool Updating           ************************/
/*****************************************************************************/
int Calculating_Distance(Solution S1, Solution S2) {
  int i;
  int dist = 0;
  for (i = 0; i < N; i++)
    if (S1.S[i] != S2.S[i]) dist += 1;
  return dist;
} // compute the distance between two solutions

int Pool_Updating_Greedy(Solution POP[], Solution S_off, int* position) {
  int i, j, k_min;
  double worst = 9999999;

  for (i = 0; i < Npop; i++)
    if (Calculating_Distance(POP[i], S_off) == 0) return 0;

  for (i = 0; i < Npop; i++)
    if (POP[i].value < worst) {
      worst = POP[i].value;
      k_min = i;
    }

  if (S_off.value > worst) {
    for (j = 0; j < N; j++)
      POP[k_min].S[j] = S_off.S[j];
    POP[k_min].value = S_off.value;
    *position = k_min;
    return 1;
  }
  return 0;
} // update the population
/*****************************************************************************/
/*****************         8. Memetic Search          ************************/
/*****************************************************************************/
void Memetic() {
  int i, j;
  int x1, x2;
  int position;
  int flag;
  int p;
  // int N_times=0;
  N_times = 0;
  SS_BEST.value = -99999999;
  starting_time = clock();
  flag = 0;

  do {

    if (flag == 0) {
      InitiaPoP(POP);
      flag = 1;
    } else
      Rebuild(POP, SS_BEST);

    time_one_run = (double)((clock() - starting_time) / CLOCKS_PER_SEC);
    if (time_one_run > time_limit) break;

    do {
      // printf("\n");
      // for(i=0;i<Npop;i++)printf("%lf  ", POP[i].value);
      // printf("\n");
      // printf("Number of pairs = %d \n",pair_s[0].number);

      p = rand() % pair_s[0].number;
      x1 = pair_s[p].i;
      x2 = pair_s[p].j;
      updating_PairSet(x1, x2);

      RandomCrossover(POP[x1], POP[x2], S_off);
      One_Flip_Tabu_Search(S_off.S, &S_off.value,
                           alpha); // printf("\n c  f = %lf \n",S_off.value);

      if (S_off.value > BestResult + 10e-7) {
        BestResult = S_off.value;
        N_times = 0;
        total_time_suc = (double)((clock() - starting_time) / CLOCKS_PER_SEC);
      } else if (fabs(S_off.value - BestResult) < 10e-7) {
        if (N_times == 0)
          total_time_suc += (double)((clock() - starting_time) / CLOCKS_PER_SEC);
        N_times++;
      }

      if (S_off.value > SS_BEST.value) {
        for (j = 0; j < N; j++)
          SS_BEST.S[j] = S_off.S[j];
        SS_BEST.value = S_off.value;
      }
      if (Pool_Updating_Greedy(POP, S_off, &position)) {
        updating_PairSet_pop(position);
      }
      time_one_run = (double)((clock() - starting_time) / CLOCKS_PER_SEC);

    } while (pair_s[0].number != 0 && time_one_run < time_limit);

    time_one_run = (double)((clock() - starting_time) / CLOCKS_PER_SEC);
  } while (time_one_run < time_limit);
}
/*****************************************************************************/
/*****************          9. Main Scheme            ************************/
/*****************************************************************************/
int main(int argc, char** argv) {
  int i, j, seed;
  int Nruns;
  // seed = time(NULL) % 1000000;
  // srand(seed);
  if (argc != 4) {
    printf("Usage: ./mammdp_fv2 {instance} {timelim} {seed}\n");
    exit(0);
  }

  File_Name = argv[1];
  time_limit = atoi(argv[2]);
  seed = atoi(argv[3]);
  // outfilename = argv[2];

  // File_Name   =  "MDPI2_500.txt";
  // outfilename = "MAMMDP_TS.txt";

  BestResult = -99999999;
  Nhit = 0;
  AvgResult = 0.0;
  total_time_suc = 0.0;

  Initializing();
  AssignMemery();

  Nruns = 1;
  // time_limit = 20;

  //  if(N<=1000)time_limit = 100;  // 100s for n<=1000, 1000s for n=3000, 2000s
  //  for n= 5000 else if(N==3000) time_limit = 1000; else time_limit = 2000;

  for (i = 0; i < Nruns; i++) {
    Memetic();
    AvgResult += SS_BEST.value;
    if (SS_BEST.value > BestResult + 10e-7) {
      BestResult = SS_BEST.value;
      Nhit = 1;
      for (j = 0; j < N; j++)
        G_BEST.S[j] = SS_BEST.S[j];
      G_BEST.value = SS_BEST.value;
    } else if (fabs(SS_BEST.value - BestResult) < 10e-7) {
      Nhit++;
      for (j = 0; j < N; j++)
        G_BEST.S[j] = SS_BEST.S[j];
      G_BEST.value = SS_BEST.value;
    }
  }
  total_time_suc /= Nhit;
  AvgResult /= Nruns;

  int bestSolSize = 0;
  for (int i = 0; i < N; ++i)
    bestSolSize += bool(G_BEST.S[i]);


  string filenameStripped(File_Name);
  // thanks https://stackoverflow.com/a/8520815
  const size_t last_slash_idx = filenameStripped.find_last_of("\\/");
  if (std::string::npos != last_slash_idx) {
    filenameStripped.erase(0, last_slash_idx + 1);
  }
  const size_t period_idx = filenameStripped.rfind('.');
  if (std::string::npos != period_idx) {
    filenameStripped.erase(period_idx);
  }

  printf("summary_line instance=%s value=%lf size=%d ttb=%lf seed=%d\n",
         filenameStripped.c_str(), BestResult, bestSolSize, total_time_suc, seed);

  // Outresulting(File_Name, outfilename, total_time_suc); // output the statistical
  // results OutSol(G_BEST, File_Name);                            // output the best
  // solution found
  return 1;
}
