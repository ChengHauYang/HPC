#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"

double GetNormal(const vector* v)
{
    const int size = v->size;
    double normsquare=0;
    for(int i=1;i<=size;i++){
      //printf("vgetp:%10.3e\n", vgetp(v,i));
      normsquare += pow(vgetp(v,i),2);
    }
    //printf("normalsquare:%10.3e\n", normsquare);
    return sqrt(normsquare);
}

vector NormalizeVecor(const vector* v)
{

  double normal = GetNormal(v);
  const int size = v->size;
  vector z = new_vector(size);

  for (int i=1; i<=size; i++)
    {
      vget(z,i) = vgetp(v,i)/normal;
    }

  return z;
}

double GetLamda(const vector* x,const matrix* A)
{
  vector Ax=matrix_vector_mult(A,x);
  double Lamda = 0;

  for (int i=1;i<=x->size;i++){
    Lamda += vget(Ax,i)*vgetp(x,i);
  }

  return Lamda;
}

double PowIter(vector* v, double TOL, int MaxIters, matrix* A)
{
    vector v_norm = NormalizeVecor(v);
    //printf("normal:%10.3e\n", normal);
    //print_vector(&v_norm);
    //print_matrix(A);

    double lambda = GetLamda(&v_norm,A);
    //printf("Lamda:%10.8e\n", lambda);
    double lambda_old = 0;
    int mstop = 0;
    int k = 0;
    while (mstop ==0){
      k++;
      vector w = matrix_vector_mult(A,&v_norm);
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
