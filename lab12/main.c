#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include <time.h> // time library

int main()
{
/*
  // Matrices
  matrix A = new_matrix(5,5);
  matrix B = new_matrix(5,5);

  for(int i=1; i<=5; i++ )
    for (int j=1; j<=5; j++ )
      {
        mget(A,i,j) = -1.0*(i==j)
                    +  2.0*(i-1==j) + 2.0*(j-1==i);
        mget(B,i,j) =  2.0*(i==j)
                    +  1.0*(i-1==j) + 1.0*(j-1==i);
      }

  // Print matrices
  print_matrix(&A);
  print_matrix(&B);

  // Add/Subtract/Multiply matrices
  matrix  Csum = matrix_add(&A,&B); print_matrix(&Csum);
  matrix Cdiff = matrix_sub(&A,&B); print_matrix(&Cdiff);
  matrix Cprod = matrix_mult(&A,&B); print_matrix(&Cprod);
  matrix  Cdot = matrix_dot_mult(&A,&B); print_matrix(&Cdot);

  // Vectors
  vector x = new_vector(5);
  vector y = new_vector(5);

  vget(x,1) = 1.0;  vget(y,1) = 1.0;
  vget(x,2) = 0.0;  vget(y,2) = 2.0;
  vget(x,3) = 1.0;  vget(y,3) = 3.0;
  vget(x,4) = 0.0;  vget(y,4) = 4.0;
  vget(x,5) = 1.0;  vget(y,5) = 5.0;

  // Print vectors
  print_vector(&x);
  print_vector(&y);

  // Add/Subtract/Multiply vectors
  vector  zsum = vector_add(&x,&y); print_vector(&zsum);
  vector zdiff = vector_sub(&x,&y); print_vector(&zdiff);
  double  zdot = vector_dot_mult(&x,&y); print_scalar(&zdot);

  // Matrix vector multiply
  vector Ax = matrix_vector_mult(&A,&x);
  print_vector(&Ax);

  // Linear solve via Gaussian elimination
  vector soln = solve(&A,&y);
  print_vector(&soln);

  // Cleanup
  delete_matrix(&A);     delete_matrix(&B);
  delete_matrix(&Csum);  delete_matrix(&Cdiff);
  delete_matrix(&Cprod); delete_matrix(&Cdot);
  delete_vector(&x);     delete_vector(&y);
  delete_vector(&zsum);  delete_vector(&zdiff);
  delete_vector(&Ax);    delete_vector(&soln);


  printf("=======================================\n");
*/

double TOL = 1e-8;
int MaxIters = 500;
double PowIter(vector* v, double TOL, int MaxIters, matrix* A);
double RQiter(vector* v, double TOL, int MaxIters, matrix* A);
vector PowIterGuessVector(vector* v, double TOL, int MaxIters, matrix* A);
vector InversePowIterGuessVector(double mu, vector* v, double TOL, int MaxIters, matrix* A);

  matrix A = new_matrix(3,3);
  vector v = new_vector(3);

  for(int i=1; i<=3; i++ ){
    mget(A,i,i) = i+1.0;
    for (int j=1; j<=3; j++ )
      {
        if (j != i){
          mget(A,i,j) = 1.0;
        }
      }
  }
  for(int i=1; i<=3; i++ ){
    vget(v,i)=1;
  }


  print_matrix(&A);
  print_vector(&v);


  printf("========Power Iteration Algorithm======== \n");
  double lambda = PowIter(&v,  TOL,  MaxIters, &A);
  printf("largest eigenvalue: %10.8e\n", lambda);

  printf("========Rayleigh Quotient Iteration Method======== \n");
  double lambda2 = RQiter(&v,  TOL,  MaxIters, &A);
  printf("largest eigenvalue: %10.8e\n", lambda2);

  printf("=======================User Input N========================== \n");


  //int N=10;

  int N;
  printf("\n Input N of matrix:");
  scanf("%i", &N);

  matrix OneTwoOne(const int N);
  matrix A1 = OneTwoOne(N);
  vector v1 = new_vector(N);
  for(int i=1; i<=N; i++ ){
    //vget(v1,i)=1; // original one
    vget(v1,i)=i; // make it not as 1
  }
  /*
  print_matrix(&A1);
  print_vector(&v1);
  */

  vector GuessEigenVector_max = PowIterGuessVector(&v1,  TOL,  N*10, &A1);
  double lambdaN_max = RQiter(&GuessEigenVector_max,  TOL,  MaxIters*N/2, &A1);
  printf("largest eigenvalue (Guess_RQ): %10.8e\n", lambdaN_max);

  //print_matrix(&A1);
  //print_vector(&v1);

  vector GuessEigenVector_min = InversePowIterGuessVector(0,&v1,  TOL,  N*10, &A1);
  double lambdaN_min = RQiter(&GuessEigenVector_min,  TOL,  MaxIters*N/2, &A1);
  printf("smallest eigenvalue (Guess_RQ): %10.8e\n", lambdaN_min);

  //print_vector(&GuessEigenVector_max);
  //print_vector(&GuessEigenVector_min);

/*
  print_matrix(&A1);
  print_vector(&v1);
*/

/*
  double lambdaN1 = RQiter(&v1,  TOL,  MaxIters, &A1);
  printf("largest eigenvalue (RQ): %10.8e\n", lambdaN1);
  double lambdaN2 = PowIter(&v1,  TOL,  MaxIters, &A1);
  printf("largest eigenvalue (POW): %10.8e\n", lambdaN2);
*/

  printf("========================Set N and calculate time========================= \n");
  double times[5];
  int k=0;

  int MatrixSize[5] = {10, 20, 40, 80, 160};

  for (int N=10;N<=160;N=N*2){
    printf("===================<N: %i>=====================\n", N);
    matrix OneTwoOne(const int N);
    matrix A1 = OneTwoOne(N);
    vector v1 = new_vector(N);
    for(int i=1; i<=N; i++ ){
      vget(v1,i)=i;
    }

/*
    print_matrix(&A1);
    print_vector(&v1);
*/

    clock_t time0, time1;
    time0 = clock();
    vector GuessEigenVector_max = PowIterGuessVector(&v1,  TOL,  N*10, &A1);
    double lambdaN_max = RQiter(&GuessEigenVector_max,  TOL,  MaxIters*N/10, &A1);
    //vector GuessEigenVector_max = PowIterGuessVector(&v1,  TOL,  50, &A1);
    //double lambdaN_max = RQiter(&GuessEigenVector_max,  TOL,  500, &A1);
    printf("largest eigenvalue (Guess_RQ): %10.8e\n", lambdaN_max);

    //print_matrix(&A1);
    vector GuessEigenVector_min = InversePowIterGuessVector(0,&v1,  TOL,  N*10, &A1);
    double lambdaN_min = RQiter(&GuessEigenVector_min,  TOL,  MaxIters*N/10, &A1);
    //vector GuessEigenVector_min = InversePowIterGuessVector(0,&v1,  TOL,  10, &A1);
    //double lambdaN_min = RQiter(&GuessEigenVector_min,  TOL,  500, &A1);
    printf("smallest eigenvalue (Guess_RQ): %10.8e\n", lambdaN_min);
    time1 = clock();

    double total_time=((double)((time1-time0))/((double)(CLOCKS_PER_SEC)));
    times[k] = total_time;
    k++;

/*
    double lambdaN1 = RQiter(&v1,  TOL,  MaxIters*N/5, &A1);
    printf("largest eigenvalue (RQ): %10.8e\n", lambdaN1);
    double lambdaN2 = PowIter(&v1,  TOL,  MaxIters*N/5, &A1);
    printf("largest eigenvalue (POW): %10.8e\n", lambdaN2);
*/
    printf("\nTotal cpu time = %10.8e s\n", total_time);

  }

  // print solution to file
  char filename[] = "cputime.data";
  FILE* outfile = fopen(filename, "w");
  for (int k=0; k<5; k++)
  {
      fprintf(outfile, "%d %11.5e\n", MatrixSize[k], times[k]);
  }
  fclose(outfile);

  // Call python script to plot
  system("python phase_plot.py");


}
