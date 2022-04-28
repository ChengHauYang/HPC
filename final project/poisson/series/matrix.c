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

matrix new_matrix(const int rows, const int cols)
{
  matrix mat;
  mat.rows = rows;
  mat.cols = cols;
  assert(rows>0);
  assert(cols>0);
  mat.val = (double*)malloc(sizeof(double)*rows*cols);

  for (int i=0; i<(rows*cols); i++)
    {  mat.val[i] = 0.0;  }

  return mat;
}

void print_matrix_full(const matrix* mat, char* varname)
{
  assert(mat->rows>0); assert(mat->cols>0);
  printf("\n %.100s =\n", &varname[1] );
  for(int i=1; i<=mat->rows; i++ )
    {
      printf("  |  ");
      for(int j=1; j<=mat->cols; j++)
        {
          printf("%10.3e",mgetp(mat,i,j));
          if (j<mat->cols) {printf(", ");}
          else {printf(" ");}
        }
      printf("|\n");
    }
  printf("\n");
}

matrix matrix_add(const matrix* A, const matrix* B)
{
  const int rows = A->rows;
  const int cols = A->cols;
  assert(rows==B->rows);
  assert(cols==B->cols);
  matrix C = new_matrix(rows,cols);

  for (int i=1; i<=rows; i++)
    for (int j=1; j<=cols; j++)
      {
        mget(C,i,j) = mgetp(A,i,j)+mgetp(B,i,j);
      }

  return C;
}

matrix matrix_sub(const matrix* A, const matrix* B)
{
  const int rows = A->rows;
  const int cols = A->cols;
  assert(rows==B->rows);
  assert(cols==B->cols);
  matrix C = new_matrix(rows,cols);

  for (int i=1; i<=rows; i++)
    for (int j=1; j<=cols; j++)
      {
        mget(C,i,j) = mgetp(A,i,j)-mgetp(B,i,j);
      }

  return C;
}

matrix scale_matrix(const matrix* A, double scaling)
{
  const int rowsA = A->rows; const int colsA = A->cols;
  matrix C = new_matrix(rowsA,colsA);

  for (int i=1; i<=colsA; i++)
      for (int k=1; k<=rowsA; k++)
        {
          mget(C,i,k) += mgetp(A,i,k)*scaling;
        }

  return C;
}

matrix matrix_mult_transpose(const matrix* A, const matrix* B)
{
  const int rowsA = A->rows; const int colsA = A->cols;
  const int colsB = B->cols;
  assert(colsA==colsB);
  matrix C = new_matrix(colsA,colsB);

  for (int i=1; i<=colsA; i++)
    for (int j=1; j<=colsB; j++)
      for (int k=1; k<=rowsA; k++)
        {
          mget(C,i,j) += mgetp(A,k,i)*mgetp(B,k,j);
        }

  return C;
}

matrix matrix_mult(const matrix* A, const matrix* B)
{
  const int rowsA = A->rows; const int colsA = A->cols;
  const int rowsB = B->rows; const int colsB = B->cols;
  assert(colsA==rowsB);
  matrix C = new_matrix(rowsA,colsB);

  for (int i=1; i<=rowsA; i++)
    for (int j=1; j<=colsB; j++)
      for (int k=1; k<=colsA; k++)
        {
          mget(C,i,j) += mgetp(A,i,k)*mgetp(B,k,j);
        }

  return C;
}

matrix matrix_dot_mult(const matrix* A, const matrix* B)
{
  const int rows = A->rows;
  const int cols = A->cols;
  assert(rows==B->rows);
  assert(cols==B->cols);
  matrix C = new_matrix(rows,cols);

  for (int i=1; i<=rows; i++)
    for (int j=1; j<=cols; j++)
      {
        mget(C,i,j) = mgetp(A,i,j)*mgetp(B,i,j);
      }

  return C;
}

void delete_matrix(matrix* mat)
{
  free(mat->val);
}

vector new_vector(const int size)
{
  vector vec;
  vec.size = size;
  assert(size>0);
  vec.val = (double*)malloc(sizeof(double)*size);

  for (int i=0; i<(size); i++)
    {  vec.val[i] = 0.0;  }

  return vec;
}

// here we return address
vector* scale_vector(const vector* x, double scaling)
{
  const int size = x->size;
  vector z = new_vector(size);

  for (int i=1; i<=size; i++)
  {
    vget(z,i) += vgetp(x,i)*scaling;
  }

  return &z;
}

vector scale_vectorTrue(const vector* x, double scaling)
{
  const int size = x->size;
  vector z = new_vector(size);

  for (int i=1; i<=size; i++)
  {
    vget(z,i) += vgetp(x,i)*scaling;
  }

  return z;
}

void print_vector_full(const vector* vec, char* varname)
{
  assert(vec->size>0);
  printf("\n");
  printf(" %.100s =\n", &varname[1] );
  printf("  |  ");
  for(int i=1; i<=vec->size; i++ )
    {
      printf("%10.3e",vgetp(vec,i));
      if (i<vec->size) {printf(", ");}
    }
  printf(" |^T\n\n");
}

vector vector_add(const vector* x, const vector* y)
{
  const int size = x->size;
  assert(size==y->size);
  vector z = new_vector(size);

  for (int i=1; i<=size; i++)
    {
      vget(z,i) = vgetp(x,i)+vgetp(y,i);
    }

  return z;
}

vector vector_sub(const vector* x, const vector* y)
{
  const int size = x->size;
  assert(size==y->size);
  vector z = new_vector(size);

  for (int i=1; i<=size; i++)
    {
      vget(z,i) = vgetp(x,i)-vgetp(y,i);
    }

  return z;
}

double vector_dot_mult(const vector* x, const vector* y)
{
  const int size = x->size; assert(size==y->size);

  double z = 0.0;
  for (int i=1; i<=size; i++)
    { z += vgetp(x,i)*vgetp(y,i); }

  return z;
}

void delete_vector(vector* vec)
{
  free(vec->val);
}

void print_scalar_full(const double* z, char* varname)
{
  printf("\n %.100s =\n", &varname[1] );
  printf("    %10.3e \n\n",*z);
}

vector matrix_vector_mult(const matrix* A, const vector* x)
{
  const int rows = A->rows; const int cols = A->cols;
  const int size = x->size;
  assert(cols==size);
  vector Ax = new_vector(rows);

  for (int i=1; i<=rows; i++)
    {
      double tmp = 0.0;
      for (int j=1; j<=size; j++)
        { tmp += mgetp(A,i,j)*vgetp(x,j); }
      vget(Ax,i) = tmp;
    }

  return Ax;
}


double det_33(const matrix* A)
{
  double r1,r2,r3,detval;

    r1 = mgetp(A,1,1) * ((mgetp(A,2,2) * mgetp(A,3,3))
    - (mgetp(A,3,2) * mgetp(A,2,3)));

    r2 = mgetp(A,1,2) * ((mgetp(A,2,1) * mgetp(A,3,3))
    - (mgetp(A,3,1) * mgetp(A,2,3)));

    r3 = mgetp(A,1,3) * ((mgetp(A,2,1) * mgetp(A,3,2))
    - (mgetp(A,3,1) * mgetp(A,2,2)));

    detval= r1 - r2 + r3;

  return detval;
}

vector solve(const matrix* A, const vector* b)
{
  const int rows = A->rows; const int cols = A->cols;
  const int size = b->size;
  assert(rows==cols); assert(rows==size);

  vector x = new_vector(rows);

  for (int i=1; i<=(size-1); i++) // LOOP OVER EACH COLUMN
    {
      // Select largest pivot in current column
      int p=i; double maxA = -100.0e0;
      for (int j=i; j<=size; j++)
        {
          double tmp = fabs(mgetp(A,j,i));
          if ( tmp > maxA ){ p = j; maxA = tmp; }
        }

      // See if matrix is singular
      if (maxA <= 1.0e-14)
        { printf(" Cannot invert system\n"); exit(1); }

      // Pivot (aka interchange rows)
      if (p!=i)
        {
          for (int j=1; j<=size; j++)
            {
              double tmp1  = mgetp(A,i,j);
              mgetp(A,i,j) = mgetp(A,p,j);
              mgetp(A,p,j) = tmp1;
            }

          double tmp2 = vgetp(b,i);
          vgetp(b,i) = vgetp(b,p);
          vgetp(b,p) = tmp2;
        }

      // Eliminate below diagonal
      for (int j=(i+1); j<=size; j++)
        {
          double dm = mgetp(A,j,i)/mgetp(A,i,i);
          for (int k=1; k<=size; k++)
            { mgetp(A,j,k) = mgetp(A,j,k) - dm*mgetp(A,i,k); }
          vgetp(b,j) = vgetp(b,j) - dm*vgetp(b,i);
        }
    }

  // Backward substitution
  vget(x,size) = vgetp(b,size)/mgetp(A,size,size);
  for (int j=1; j<=(size-1); j++)
    {
      double sum = 0.0e0;

      for (int k=(size-j+1); k<=size; k++)
        { sum = sum + mgetp(A,size-j,k)*vget(x,k); }

      vget(x,size-j) = (vgetp(b,size-j) - sum)
                     /mgetp(A,size-j,size-j);
    }

  return x;
}


vector solveCG(const matrix* A, const vector* b)
{
  const int rows = A->rows; const int cols = A->cols;
  const int size = b->size;
  assert(rows==cols); assert(rows==size);

  vector x = new_vector(size);
  vector r = *b;
  double rho = vector_dot_mult(&r,&r);

  printf("rho = %25.15f \n",rho);

  vector p = r;
  double residual = GetNormal(&r);

  printf("residual = %25.15f \n",residual);

  vector q = new_vector(size);
  vector temp = new_vector(size);  
  double alpha,rho_old, beta,residual_old,residual_ori;

  residual_ori = residual;
  for (int iter=1;iter<=5000;iter++){
    q = matrix_vector_mult(A,&p);
    alpha = rho/(vector_dot_mult(&p,&q));
    printf("alpha = %25.15f \n",alpha);

    //x=vector_add(&x,scale_vector(&p,alpha));
    //r=vector_add(&r,scale_vector(&q,-alpha));
    temp=scale_vectorTrue(&p,alpha);
    x=vector_add(&x,&temp);
    temp=scale_vectorTrue(&q,-alpha);
    r=vector_add(&r,&temp);    
    rho_old = rho;
    rho =vector_dot_mult(&r,&r);
    beta = rho / rho_old;
    printf("beta = %25.15f \n",beta);

    //p=vector_add(&r,scale_vector(&p,beta));
    temp=scale_vectorTrue(&p,beta);
    p=vector_add(&r,&temp);

    residual_old = residual;
    residual = GetNormal(&r);
    if (residual/residual_ori < 1e-13){
      break;
    }
  }

  return x;
}

/*
vector solveCGpre(const matrix* A, const vector* b)
{
  const int rows = A->rows; const int cols = A->cols;
  const int size = b->size;
  assert(rows==cols); assert(rows==size);

  vector x = new_vector(size);
  vector r = *b;
  vector z = r; // need to change if we have preconditioner
  double rho = vector_dot_mult(&r,&z);
  vector p = z;
  double residual = GetNormal(&r);


  vector q;
  double alpha,rho_old, beta,residual_old,residual_ori;

  residual_ori = residual;
  for (int iter=1;iter<=5000;iter++){
    q = matrix_vector_mult(A,&p);
    alpha = rho/(vector_dot_mult(&p,&q));
    x=vector_add(&x,scale_vector(&p,alpha));
    r=vector_add(&r,scale_vector(&q,-alpha));
    z = r;  // need to change if we have preconditioner
    rho_old = rho;
    rho =vector_dot_mult(&r,&z);
    beta = rho / rho_old;
    p=vector_add(&z,scale_vector(&p,beta));

    residual_old = residual;
    residual = GetNormal(&r);
    if (residual/residual_ori < 1e-13){
      break;
    }
  }

  return x;
}
*/
