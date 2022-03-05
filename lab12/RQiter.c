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

matrix ObtainLHS(const matrix* A,const int lamda)
{
  const int col = A->cols;
  matrix C = new_matrix(col,col);

  for (int i=1; i<=col; i++){
    for (int j=1; j<=col; j++)
      {
        if (j != i){
          mget(C,i,j) = mgetp(A,i,j);
        } else
          mget(C,i,i) = mgetp(A,i,j)-lamda;
      }
  }

  return C;
}


double RQiter(vector* v, double TOL, int MaxIters, matrix* A)
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
      matrix LHS = ObtainLHS(A,lambda);
      vector w = solve(&LHS,&v_norm);
      v_norm = NormalizeVecor(&w);
      if ((fabs(lambda-lambda_old) < TOL) || (k==MaxIters)){
        mstop =1;
      }
      lambda_old = lambda;
      lambda = GetLamda(&v_norm,A);
      //printf("%10.8e\n", lambda);
    }

    printf("Iterations: %d\n", k);
    //printf("Lamda:%10.8e\n", lambda);

    return lambda;
}
