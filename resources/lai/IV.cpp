  
/*****************************************************************************\
 * This program aims to generate the benchmark instances used in our paper   *
 * for the weighted Max-Mean dispersion problem.                             * 
\*****************************************************************************/ 
                                                                                                                                         
/*****************************************************************************/
/**********          0. Header files and variables      **********************/
/*****************************************************************************/ 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <math.h>
// #include <conio.h>
#include <ctype.h>
using namespace std;
FILE * fp ;
int N;        // node number in graph
double  ** D; // distance matrix
double  * w ; // node weight 
const double MW = 1.0; 
/*****************************************************************************/
/*****************          1. Outputing  results      ***********************/
/*****************************************************************************/ 
void Outputing(int tt)
{
    int i,j;
	FILE *fp; 
	char buff[80];
	
    sprintf(buff,"%s_%d_%d.txt","IV",N,tt);
    fp=fopen(buff,"a+");
    fprintf(fp,"%d\n", N); 
    for(i=0;i<N;i++)  fprintf(fp,"%d   %5.2lf\n",i+1, w[i]); 
    for(i=0;i<N;i++)
    for(j=i+1;j<N;j++) 
      fprintf(fp,"%d   %d   %5.2lf\n",i+1, j+1, D[i][j]); 
	fclose(fp);
}

/*****************************************************************************/
/*****************         2. Main Function           ************************/
/*****************************************************************************/ 
int main(int argc, char **argv)
{ 
     int i, j, i1, j1, seed ; 
     int tt;

     N = atoi(argv[1]);
     seed = N;

    //  seed = 5000; //seed= 1000 for N = 1000, seed= 3000 for N =3000, seed= 5000 for N = 5000
    //  srand( seed ) ;
  
     int sign;
     N = 5000;    //N = 1000, 3000, 5000. 
     w = new double [N]; 
     D = new double *[N];
     for(i=0;i<N;i++)
     D[i] = new double [N];
     for(tt=1;tt<=10;tt++)
     {
      for(i=0;i<N;i++)
      for(j=i+1;j<N;j++)
      {  
         D[i][j] = (-1 + rand()%3)*10.0 ; 
         D[j][i] = D[i][j]; 
      } 
      
      for(i=0;i<N;i++) 
      {  
        w[i] = 1.0 ;
      }
      
      Outputing(tt);
    }

  return 0;
}
 
