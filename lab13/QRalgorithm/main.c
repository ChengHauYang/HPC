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
    if (T->center[i] < *largest_abs_eig){
      *smallest_abs_eig=T->center[i];
    }
  }

  return *largest_abs_eig/ *smallest_abs_eig;
}



matrix OneTwoOne(const int N)
{
  matrix A = new_matrix(N,N);
  double adding = 1/(double)N/(double)N;
  //printf("%d\n", N);
  //printf("%10.8e \n", adding);

  for(int i=1; i<=N; i++ )
    for (int j=1; j<=N; j++ )
      {
        mget(A,i,j) = (2.0+adding)*(i==j)
                    -  1.0*(i-1==j) - 1.0*(j-1==i);
      }
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

  // Hessenberg reduction to tridiagonal form
  trimatrix T = new_trimatrix(N);
  void Hessenberg(const matrix* A, trimatrix* T);
  Hessenberg(&A,&T);
  printf("\n");
  printf("Original Matrix:\n");
  print_matrix(&A);
  printf("Reduction to Tridiagonal Form:\n");
  print_trimatrix(&T);

  // QR Algorithm to find eigenvalues of T
  // which are the same as the eigenvalues of A
  void QRA(trimatrix* T);
  QRA(&T);
  printf("After QR Algorithm:\n");
  print_trimatrix(&T);

  double smallest_abs_eig,largest_abs_eig;

  double conditionnumber=ConditionNumber(&T,&smallest_abs_eig,&largest_abs_eig);

  printf("largest eigenvalue:%10.3e\n",largest_abs_eig);
  printf("smallest eigenvalue:%10.3e\n",smallest_abs_eig);
  printf("Condition Number:%10.3e\n",conditionnumber);



}
