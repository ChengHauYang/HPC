#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "matrix.h"
#include "trimatrix.h"

double ConditionNumber(const trimatrix* T,
                       double* smallest_abs_eig,
                       double* largest_abs_eig)
{
/*
Print the smallest eigenvalue, largest eigenvalue, and the condition number to the screen.
*/

*smallest_abs_eig= T->center[0];
*largest_abs_eig= T->center[0];
for(int i=1; i<T->rows; i++ )
  {
    if (T->center[i] > *largest_abs_eig){
      *largest_abs_eig=T->center[i];
    }
    if (T->center[i] < *smallest_abs_eig){
      *smallest_abs_eig=T->center[i];
    }
  }

  return *largest_abs_eig/ *smallest_abs_eig;
}



matrix OneTwoOne(const int N)
{
  matrix A = new_matrix(N,N);
  double adding = 1/(double)N/(double)N;
  double mult = (double)N*(double)N;
  //printf("%d\n", N);
  //printf("%10.8e \n", adding);

  for(int i=1; i<=N; i++ )
    for (int j=1; j<=N; j++ )
      {
        mget(A,i,j) = mult*(2.0+adding)*(i==j)
                    -  mult*1.0*(i-1==j) - mult*1.0*(j-1==i);
      }
  mget(A,N,1) = - mult*1.0;
  mget(A,1,N) = - mult*1.0;
  return A;
}

int main()
{
  srand( time(NULL) );
  int N=0;
  printf("    Input N: ");
  scanf("%i",&N);
  assert(N>1);

  // Create a matrix
  matrix A = OneTwoOne(N);

/*
  for (int i=1; i<=N; i++)
    for (int j=1; j<=i; j++)
      {
        double tmp = ((double)rand())/((double)RAND_MAX);
        tmp = (double)((int)(10000.0*tmp))/10000.0;
        mget(A,i,j) = tmp;
      }
  for (int i=1; i<=(N-1); i++)
    for (int j=(i+1); j<=N; j++)
      { mget(A,i,j) = mget(A,j,i); }
*/

  clock_t time0, time1, time2, time3, time4, time5;

  // Hessenberg reduction to tridiagonal form
  trimatrix T = new_trimatrix(N);
  void Hessenberg(const matrix* A, trimatrix* T);

  time0 = clock();
  Hessenberg(&A,&T);
  time1 = clock();

  printf("\n");
  printf("Original Matrix:\n");
  print_matrix(&A);
  printf("Reduction to Tridiagonal Form:\n");
  print_trimatrix(&T);

  // QR Algorithm to find eigenvalues of T
  // which are the same as the eigenvalues of A

  void QRA(trimatrix* T);

  time2 = clock();
  QRA(&T);
  time3 = clock();

  printf("After QR Algorithm:\n");
  print_trimatrix(&T);

  double smallest_abs_eig,largest_abs_eig;

  time4 = clock();
  double conditionnumber=ConditionNumber(&T,&smallest_abs_eig,&largest_abs_eig);
  time5 = clock();

  double TimeHessenberg=((double)((time1-time0))/((double)(CLOCKS_PER_SEC)));
  double TimeQRA=((double)((time3-time2))/((double)(CLOCKS_PER_SEC)));
  double TimeCondi=((double)((time5-time4))/((double)(CLOCKS_PER_SEC)));

  printf("largest eigenvalue:%10.3e\n",largest_abs_eig);
  printf("smallest eigenvalue:%10.3e\n",smallest_abs_eig);
  printf("Condition Number:%10.3e\n",conditionnumber);

  printf("Hessenberg time = %10.8e s\n", TimeHessenberg);
  printf("QRA time = %10.8e s\n", TimeQRA);
  printf("ConditionNumber time = %10.8e s\n", TimeCondi);


////////////////////////

printf("========================Set N and calculate time========================= \n");

int MatrixSize[5] = {10, 20, 40, 80, 160};
double TimesHessenberg[5];
double TimesQRA[5];
double TimesCondi[5];
double conditionnumbers[5];

for (int i=0;i<5;i++){
//for (int i=2;i<3;i++){
  printf("\n===================<N: %i>=====================\n", MatrixSize[i]);
  clock_t time0, time1, time2, time3, time4, time5;

  matrix A = OneTwoOne(MatrixSize[i]);

  // Hessenberg reduction to tridiagonal form
  trimatrix T = new_trimatrix(MatrixSize[i]);

  time0 = clock();
  Hessenberg(&A,&T);
  time1 = clock();

/*
  printf("\n");
  printf("Original Matrix:\n");
  print_matrix(&A);
*/

/*
  printf("\n");
  printf("Original Matrix:\n");
  print_matrix(&A);
  printf("Reduction to Tridiagonal Form:\n");
  print_trimatrix(&T);
*/

  // QR Algorithm to find eigenvalues of T
  // which are the same as the eigenvalues of A

  time2 = clock();
  QRA(&T);
  time3 = clock();

/*
  printf("After QR Algorithm:\n");
  print_trimatrix(&T);
*/

  double smallest_abs_eig,largest_abs_eig;


  time4 = clock();
  double conditionnumber=ConditionNumber(&T,&smallest_abs_eig,&largest_abs_eig);
  time5 = clock();

  double TimeHessenberg=((double)((time1-time0))/((double)(CLOCKS_PER_SEC)));
  double TimeQRA=((double)((time3-time2))/((double)(CLOCKS_PER_SEC)));
  double TimeCondi=((double)((time5-time4))/((double)(CLOCKS_PER_SEC)));

  printf("largest eigenvalue:%10.3e\n",largest_abs_eig);
  printf("smallest eigenvalue:%10.3e\n",smallest_abs_eig);
  printf("Condition Number:%10.3e\n",conditionnumber);

  printf("Hessenberg time = %15.8e s\n", TimeHessenberg);
  printf("QRA time = %15.8e s\n", TimeQRA);
  printf("ConditionNumber time = %25.16e s\n", TimeCondi);

  TimesHessenberg[i]=TimeHessenberg;
  TimesQRA[i]=TimeQRA;
  TimesCondi[i]=TimeCondi;
  conditionnumbers[i]=conditionnumber;
}

  printf("========================Scaling with N========================= \n");
  printf("Hessenberg time scale with N:%10.3e\n",(log(TimesHessenberg[4])-log(TimesHessenberg[0]))/(log(MatrixSize[4])-log(MatrixSize[0])));
  printf("QRA time scale with N:%10.3e\n",(log(TimesQRA[4])-log(TimesQRA[0]))/(log(MatrixSize[4])-log(MatrixSize[0])));
  printf("ConditionNumber time scale with N:%10.3e\n",(log(TimesCondi[2])-log(TimesCondi[0]))/(log(MatrixSize[2])-log(MatrixSize[0])));
  printf("conditionnumber scale with N:%10.3e\n",(log(conditionnumbers[4])-log(conditionnumbers[0]))/(log(MatrixSize[4])-log(MatrixSize[0])));


  // print solution to file
  char filename[] = "cputime.data";
  FILE* outfile = fopen(filename, "w");
  for (int k=0; k<5; k++)
  {
      fprintf(outfile, "%d %15.8e %15.8e %25.16e %15.8e\n", MatrixSize[k], TimesHessenberg[k], TimesQRA[k], TimesCondi[k],conditionnumbers[k]);
  }
  fclose(outfile);



}
