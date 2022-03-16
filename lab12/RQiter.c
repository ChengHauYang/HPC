#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"


matrix IdentityMatrix(const int num, const int col)
{
  matrix C = new_matrix(col,col);

  for (int j=1; j<=col; j++)
  {
    mget(C,j,j) = num;
  }

  return C;
}

matrix ObtainLHS(const matrix* A,const double lamda)
{
  const int col = A->cols;
  matrix C = new_matrix(col,col);

  for (int i=1; i<=col; i++){
    for (int j=1; j<=col; j++)
      {
        if (j != i){
          mget(C,i,j) = mgetp(A,i,j);
        } else
          mget(C,i,i) = mgetp(A,i,i)-lamda;
      }
  }

  return C;
}

matrix OneTwoOne(const int N)
{
  matrix A = new_matrix(N,N);
  for(int i=1; i<=N; i++ )
    for (int j=1; j<=N; j++ )
      {
        mget(A,i,j) = 2.0*(i==j)
                    -  1.0*(i-1==j) - 1.0*(j-1==i);
      }
  return A;
}

double RQiter(const vector* v, double TOL, int MaxIters,const matrix* A)
{
    vector NormalizeVecor(const vector* v);
    vector v_norm = NormalizeVecor(v);
    //printf("normal:%10.3e\n", normal);
    //print_vector(&v_norm);
    //print_matrix(A);

    double GetLamda(const vector* x,const matrix* A);
    double lambda = GetLamda(&v_norm,A);
    //printf("Lamda:%10.8e\n", lambda);
    double lambda_old = 0;
    int mstop = 0;
    int k = 0;

    //matrix I_test=IdentityMatrix(3,3);
    //print_matrix(&I_test);

    while (mstop ==0){
      k++;
      //printf("lambda:%10.3e\n", lambda);
      matrix LHS = ObtainLHS(A,lambda);
      //print_matrix(&LHS);
      vector w = solve(&LHS,&v_norm);
      v_norm = NormalizeVecor(&w);
      //print_vector(&v_norm);
      lambda_old = lambda;
      lambda = GetLamda(&v_norm,A);
      if ((fabs(lambda-lambda_old) < TOL) || (k==MaxIters)){
        mstop =1;
      }
      //printf("%10.8e\n", lambda);
    }

    printf("Iterations: %d\n", k);
    //printf("Lamda:%10.8e\n", lambda);

    return lambda;
}
