#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include <mpi.h>
#include <vector>
#include <iostream>

matrix new_matrix(const int rows, const int cols)
{
  matrix mat;
  mat.rows = rows;
  mat.cols = cols;
  assert(rows > 0);
  assert(cols > 0);
  mat.val = (double *)malloc(sizeof(double) * rows * cols);

  for (int i = 0; i < (rows * cols); i++)
  {
    mat.val[i] = 0.0;
  }

  return mat;
}

void print_matrix_full(const matrix *mat, char *varname)
{
  assert(mat->rows > 0);
  assert(mat->cols > 0);
  printf("\n %.100s =\n", &varname[1]);
  for (int i = 1; i <= mat->rows; i++)
  {
    printf("  |  ");
    for (int j = 1; j <= mat->cols; j++)
    {
      printf("%10.3e", mgetp(mat, i, j));
      if (j < mat->cols)
      {
        printf(", ");
      }
      else
      {
        printf(" ");
      }
    }
    printf("|\n");
  }
  printf("\n");
}

matrix matrix_add(const matrix *A, const matrix *B)
{
  const int rows = A->rows;
  const int cols = A->cols;
  assert(rows == B->rows);
  assert(cols == B->cols);
  matrix C = new_matrix(rows, cols);

  for (int i = 1; i <= rows; i++)
    for (int j = 1; j <= cols; j++)
    {
      mget(C, i, j) = mgetp(A, i, j) + mgetp(B, i, j);
    }

  return C;
}

matrix matrix_sub(const matrix *A, const matrix *B)
{
  const int rows = A->rows;
  const int cols = A->cols;
  assert(rows == B->rows);
  assert(cols == B->cols);
  matrix C = new_matrix(rows, cols);

  for (int i = 1; i <= rows; i++)
    for (int j = 1; j <= cols; j++)
    {
      mget(C, i, j) = mgetp(A, i, j) - mgetp(B, i, j);
    }

  return C;
}

matrix scale_matrix(const matrix *A, double scaling)
{
  const int rowsA = A->rows;
  const int colsA = A->cols;
  matrix C = new_matrix(rowsA, colsA);

  for (int i = 1; i <= colsA; i++)
    for (int k = 1; k <= rowsA; k++)
    {
      mget(C, i, k) += mgetp(A, i, k) * scaling;
    }

  return C;
}

int find_min_matrix(const matrix *A)
{
  const int rowsA = A->rows;
  const int colsA = A->cols;
  matrix C = new_matrix(rowsA, colsA);

  int min = mgetp(A, 1, 1);
  for (int i = 1; i <= colsA; i++)
    for (int k = 1; k <= rowsA; k++)
    {
      if (mgetp(A, i, k) < min)
      {
        min = mgetp(A, i, k);
      }
    }

  return min;
}

matrix matrix_mult_transpose(const matrix *A, const matrix *B)
{
  const int rowsA = A->rows;
  const int colsA = A->cols;
  const int colsB = B->cols;
  assert(colsA == colsB);
  matrix C = new_matrix(colsA, colsB);

  for (int i = 1; i <= colsA; i++)
    for (int j = 1; j <= colsB; j++)
      for (int k = 1; k <= rowsA; k++)
      {
        mget(C, i, j) += mgetp(A, k, i) * mgetp(B, k, j);
      }

  return C;
}

matrix matrix_mult(const matrix *A, const matrix *B)
{
  const int rowsA = A->rows;
  const int colsA = A->cols;
  const int rowsB = B->rows;
  const int colsB = B->cols;
  assert(colsA == rowsB);
  matrix C = new_matrix(rowsA, colsB);

  for (int i = 1; i <= rowsA; i++)
    for (int j = 1; j <= colsB; j++)
      for (int k = 1; k <= colsA; k++)
      {
        mget(C, i, j) += mgetp(A, i, k) * mgetp(B, k, j);
      }

  return C;
}

matrix matrix_dot_mult(const matrix *A, const matrix *B)
{
  const int rows = A->rows;
  const int cols = A->cols;
  assert(rows == B->rows);
  assert(cols == B->cols);
  matrix C = new_matrix(rows, cols);

  for (int i = 1; i <= rows; i++)
    for (int j = 1; j <= cols; j++)
    {
      mget(C, i, j) = mgetp(A, i, j) * mgetp(B, i, j);
    }

  return C;
}

void delete_matrix(matrix *mat)
{
  free(mat->val);
}

vector new_vector(const int size)
{
  vector vec;
  vec.size = size;
  assert(size > 0);
  vec.val = (double *)malloc(sizeof(double) * size);

  for (int i = 0; i < (size); i++)
  {
    vec.val[i] = 0.0;
  }

  return vec;
}

// here we return address
vector *scale_vector(const vector *x, double scaling)
{
  const int size = x->size;
  vector z = new_vector(size);

  for (int i = 1; i <= size; i++)
  {
    vget(z, i) += vgetp(x, i) * scaling;
  }

  return &z;
}

vector scale_vectorTrue(const vector *x, double scaling)
{
  const int size = x->size;
  vector z = new_vector(size);

  for (int i = 1; i <= size; i++)
  {
    vget(z, i) += vgetp(x, i) * scaling;
  }

  return z;
}

void print_vector_full(const vector *vec, char *varname)
{
  assert(vec->size > 0);
  printf("\n");
  printf(" %.100s =\n", &varname[1]);
  printf("  |  ");
  for (int i = 1; i <= vec->size; i++)
  {
    printf("%10.3e", vgetp(vec, i));
    if (i < vec->size)
    {
      printf(", ");
    }
  }
  printf(" |^T\n\n");
}

vector vector_add(const vector *x, const vector *y)
{
  const int size = x->size;
  assert(size == y->size);
  vector z = new_vector(size);

  for (int i = 1; i <= size; i++)
  {
    vget(z, i) = vgetp(x, i) + vgetp(y, i);
  }

  return z;
}

vector vector_sub(const vector *x, const vector *y)
{
  const int size = x->size;
  assert(size == y->size);
  vector z = new_vector(size);

  for (int i = 1; i <= size; i++)
  {
    vget(z, i) = vgetp(x, i) - vgetp(y, i);
  }

  return z;
}

double vector_dot_mult(const vector *x, const vector *y)
{
  const int size = x->size;
  assert(size == y->size);

  double z = 0.0;
  for (int i = 1; i <= size; i++)
  {
    z += vgetp(x, i) * vgetp(y, i);
  }

  return z;
}

void delete_vector(vector *vec)
{
  free(vec->val);
}

void print_scalar_full(const double *z, char *varname)
{
  printf("\n %.100s =\n", &varname[1]);
  printf("    %10.3e \n\n", *z);
}

vector matrix_vector_mult(const matrix *A, const vector *x)
{
  const int rows = A->rows;
  const int cols = A->cols;
  const int size = x->size;
  assert(cols == size);
  vector Ax = new_vector(rows);

  for (int i = 1; i <= rows; i++)
  {
    double tmp = 0.0;
    for (int j = 1; j <= size; j++)
    {
      tmp += mgetp(A, i, j) * vgetp(x, j);
    }
    vget(Ax, i) = tmp;
  }

  return Ax;
}

double det_33(const matrix *A)
{
  double r1, r2, r3, detval;

  r1 = mgetp(A, 1, 1) * ((mgetp(A, 2, 2) * mgetp(A, 3, 3)) - (mgetp(A, 3, 2) * mgetp(A, 2, 3)));

  r2 = mgetp(A, 1, 2) * ((mgetp(A, 2, 1) * mgetp(A, 3, 3)) - (mgetp(A, 3, 1) * mgetp(A, 2, 3)));

  r3 = mgetp(A, 1, 3) * ((mgetp(A, 2, 1) * mgetp(A, 3, 2)) - (mgetp(A, 3, 1) * mgetp(A, 2, 2)));

  detval = r1 - r2 + r3;

  return detval;
}

vector solve(const matrix *A, const vector *b)
{
  const int rows = A->rows;
  const int cols = A->cols;
  const int size = b->size;
  assert(rows == cols);
  assert(rows == size);

  vector x = new_vector(rows);

  for (int i = 1; i <= (size - 1); i++) // LOOP OVER EACH COLUMN
  {
    // Select largest pivot in current column
    int p = i;
    double maxA = -100.0e0;
    for (int j = i; j <= size; j++)
    {
      double tmp = fabs(mgetp(A, j, i));
      if (tmp > maxA)
      {
        p = j;
        maxA = tmp;
      }
    }

    // See if matrix is singular
    if (maxA <= 1.0e-14)
    {
      printf(" Cannot invert system\n");
      exit(1);
    }

    // Pivot (aka interchange rows)
    if (p != i)
    {
      for (int j = 1; j <= size; j++)
      {
        double tmp1 = mgetp(A, i, j);
        mgetp(A, i, j) = mgetp(A, p, j);
        mgetp(A, p, j) = tmp1;
      }

      double tmp2 = vgetp(b, i);
      vgetp(b, i) = vgetp(b, p);
      vgetp(b, p) = tmp2;
    }

    // Eliminate below diagonal
    for (int j = (i + 1); j <= size; j++)
    {
      double dm = mgetp(A, j, i) / mgetp(A, i, i);
      for (int k = 1; k <= size; k++)
      {
        mgetp(A, j, k) = mgetp(A, j, k) - dm * mgetp(A, i, k);
      }
      vgetp(b, j) = vgetp(b, j) - dm * vgetp(b, i);
    }
  }

  // Backward substitution
  vget(x, size) = vgetp(b, size) / mgetp(A, size, size);
  for (int j = 1; j <= (size - 1); j++)
  {
    double sum = 0.0e0;

    for (int k = (size - j + 1); k <= size; k++)
    {
      sum = sum + mgetp(A, size - j, k) * vget(x, k);
    }

    vget(x, size - j) = (vgetp(b, size - j) - sum) / mgetp(A, size - j, size - j);
  }

  return x;
}

double GlobalSum(const int my_rank,
                 const int comm_sz,
                 const double local_in)
{
  double global_out;

  if (my_rank != 0)
  {
    MPI_Send(&local_in, 1, MPI_DOUBLE, 0, 200, MPI_COMM_WORLD);
  }
  else
  {
    global_out = local_in;
    for (int i = 1; i < comm_sz; i++)
    {
      double tmp_val;
      MPI_Recv(&tmp_val, 1, MPI_DOUBLE, i, 200,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      global_out += tmp_val;
    }
    for (int i = 1; i < comm_sz; i++)
    {
      MPI_Send(&global_out, 1, MPI_DOUBLE, i, 300, MPI_COMM_WORLD);
    }
  }

  if (my_rank != 0)
  {
    MPI_Recv(&global_out, 1, MPI_DOUBLE, 0, 300,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  return global_out;
}

double dot_prod(const vector *vec1, const vector *vec2)
{
  const int N = vec1->size;

  double dot = 0.0;
  for (int i = 1; i <= N; i++)
  {
    dot += vgetp(vec1, i) * vgetp(vec2, i);
  }

  return dot;
}

double stdvector_dot_mult(const std::vector<double> &x, const std::vector<double> &y)
{
  int size = x.size();
  double z = 0.0;
  for (int i = 0; i < size; i++)
  {
    z += x[i] * y[i];
  }

  return z;
}

double stdGetPower(const std::vector<double> &v)
{
  const int size = v.size();
  double normsquare = 0;
  for (int i = 0; i < size; i++)
  {
    normsquare += pow(v[i], 2);
  }
  // printf("normalsquare:%10.3e\n", normsquare);
  return (normsquare);
}

std::vector<double> scale_stdvector(const std::vector<double> &x, double scaling)
{
  const int size = x.size();
  std::vector<double> z(size);

  for (int i = 0; i < size; i++)
  {
    z[i] = x[i] * scaling;
  }

  return z;
}

std::vector<double> stdvector_add(const std::vector<double> &x, const std::vector<double> &y)
{
  const int size = x.size();
  assert(size == y.size());
  std::vector<double> z(size);

  for (int i = 0; i < size; i++)
  {
    z[i] = x[i] + y[i];
  }

  return z;
}

void send_boundary_data2D(const int my_rank, const int comm_sz,
                          const std::vector<double> &Uold)
{
  const int last_rank = comm_sz - 1;
  if (last_rank == 0)
  {
    return;
  }

  int old_size = Uold.size();

  if (my_rank == 0)
  {
    MPI_Send(&old_size, 1, MPI_INT, 1, 999, MPI_COMM_WORLD);
    MPI_Send(&Uold[0], Uold.size(), MPI_DOUBLE, 1, 999, MPI_COMM_WORLD);
  }
  else if (my_rank == last_rank)
  {
    MPI_Send(&old_size, 1, MPI_INT, last_rank - 1, 999, MPI_COMM_WORLD);
    MPI_Send(&Uold[0], Uold.size(), MPI_DOUBLE, last_rank - 1, 999,
             MPI_COMM_WORLD);
  }
  else
  {
    MPI_Send(&old_size, 1, MPI_INT, my_rank - 1, 999, MPI_COMM_WORLD);
    MPI_Send(&old_size, 1, MPI_INT, my_rank + 1, 999, MPI_COMM_WORLD);

    MPI_Send(&Uold[0], Uold.size(), MPI_DOUBLE, my_rank - 1, 999,
             MPI_COMM_WORLD);
    MPI_Send(&Uold[0], Uold.size(), MPI_DOUBLE, my_rank + 1, 999,
             MPI_COMM_WORLD);
  }
}

template <typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n)
{
  int first = v.cbegin() + m;
  int last = v.cbegin() + n + 1;

  std::vector<T> vec(first, last);
  return vec;
}

void send_boundary_data2D_new(const int my_rank, const int comm_sz,
                              const std::vector<double> &Uold, const int m)
{
  const int last_rank = comm_sz - 1;
  if (last_rank == 0)
  {
    return;
  }

  int old_size = m;

  std::vector<double> Utop(Uold.end() - m, Uold.end());
  std::vector<double> Ubot = slice(Uold, 0, m - 1);

  // std::cout<< "m=" << m <<"\n";
  // std::cout<< "utop size:" << Utop.size() <<"\n";
  // std::cout<< "ubot size:" << Ubot.size() <<"\n";

  if (my_rank == 0)
  {
    MPI_Send(&old_size, 1, MPI_INT, 1, 999, MPI_COMM_WORLD);
    MPI_Send(&Utop[0], Utop.size(), MPI_DOUBLE, 1, 999, MPI_COMM_WORLD);
  }
  else if (my_rank == last_rank)
  {
    MPI_Send(&old_size, 1, MPI_INT, last_rank - 1, 999, MPI_COMM_WORLD);
    MPI_Send(&Ubot[0], Ubot.size(), MPI_DOUBLE, last_rank - 1, 999,
             MPI_COMM_WORLD);
  }
  else
  {
    MPI_Send(&old_size, 1, MPI_INT, my_rank - 1, 999, MPI_COMM_WORLD);
    MPI_Send(&old_size, 1, MPI_INT, my_rank + 1, 999, MPI_COMM_WORLD);

    MPI_Send(&Ubot[0], Ubot.size(), MPI_DOUBLE, my_rank - 1, 999,
             MPI_COMM_WORLD);
    MPI_Send(&Utop[0], Utop.size(), MPI_DOUBLE, my_rank + 1, 999,
             MPI_COMM_WORLD);
  }
}

void receive_boundary_data2D(const int my_rank,
                             const int comm_sz,
                             std::vector<double> &Ubot,
                             std::vector<double> &Utop)
{
  const int last_rank = comm_sz - 1;
  // std::cout<<"last_rank:"<<last_rank<<"\n";
  Utop = {};
  Ubot = {};

  int top_size;
  int bot_size;
  if (last_rank == 0)
  {
    return;
  }
  if (my_rank == 0)
  {
    MPI_Recv(&top_size, 1, MPI_INT, 1, 999,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Utop.resize(top_size);
    MPI_Recv(&Utop[0], top_size, MPI_DOUBLE, 1, 999,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank == last_rank)
  {
    MPI_Recv(&bot_size, 1, MPI_INT, last_rank - 1, 999,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Ubot.resize(bot_size);
    MPI_Recv(&Ubot[0], bot_size, MPI_DOUBLE, last_rank - 1, 999,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else
  {
    MPI_Recv(&bot_size, 1, MPI_INT, my_rank - 1, 999,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&top_size, 1, MPI_INT, my_rank + 1, 999,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    Utop.resize(top_size);
    Ubot.resize(bot_size);

    MPI_Recv(&Ubot[0], bot_size, MPI_DOUBLE, my_rank - 1, 999,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&Utop[0], top_size, MPI_DOUBLE, my_rank + 1, 999,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void MatMult2D(const int my_rank, const int comm_sz,
               const std::vector<int> &Pos_local_x,
               const std::vector<int> &Pos_local_y,
               const std::vector<double> &Number_local,
               const std::vector<double> &Uold,
               std::vector<double> &U,
               const int m)
{

  void print_double_vector_seq(const int comm_sz, const int my_rank, std::vector<double> const &input);

  void send_boundary_data2D_new(const int my_rank, const int comm_sz,
                                const std::vector<double> &Uold, const int m);

  void send_boundary_data2D(const int my_rank, const int comm_sz,
                            const std::vector<double> &Uold);

  void receive_boundary_data2D(const int my_rank,
                               const int comm_sz,
                               std::vector<double> &Ubot,
                               std::vector<double> &Utop);

  const int N = Uold.size();

  // printf("Uold:");
  // print_double_vector_seq(comm_sz,my_rank,Uold);

  U.resize(N);

  std::vector<double> Utop;
  std::vector<double> Ubot;

  //  std::vector<double> Umerge;

  // send_boundary_data2D(my_rank, comm_sz, Uold);
  send_boundary_data2D_new(my_rank, comm_sz, Uold, m);
  receive_boundary_data2D(my_rank, comm_sz, Ubot, Utop);

  // only for c++ 17
  // std::merge(Ubot.begin(), Ubot.end(), Uold.begin(), Uold.end(), std::back_inserter(Umerge));
  // std::merge(Umerge.begin(), Umerge.end(), Utop.begin(), Utop.end(), std::back_inserter(Umerge));

  Ubot.insert(Ubot.end(), Uold.begin(), Uold.end());
  Ubot.insert(Ubot.end(), Utop.begin(), Utop.end());

  // printf("Ubot:");
  // print_double_vector_seq(comm_sz,my_rank,Ubot);

  for (int i = 0; i < N; i++)
  {
    // initalization
    U[i] = 0;
    for (int a = 0; a < Pos_local_x.size(); a++)
    {
      if (Pos_local_y[a] == i + 1) // Pos_local_y starts from 1 !!!!!!!!!!!!!!!
      {
        // NOTE: be careful here
        U[i] += Number_local[a] * Ubot[Pos_local_x[a] - 1]; // Pos_local_x starts from 1 !!!!!!!!!!!!!!!
      }
    }
  }

  // printf("U:");
  // print_double_vector_seq(comm_sz,my_rank,U);
}

std::vector<double> solveCGMPI(const int my_rank, const int comm_sz,
                               const std::vector<double> &F_local,
                               const std::vector<int> &Pos_local_x,
                               const std::vector<int> &Pos_local_y,
                               const std::vector<double> &Number_local,
                               const int m)
{

  double GlobalSum(const int my_rank,
                   const int comm_sz,
                   const double local_in);

  void MatMult2D(const int my_rank, const int comm_sz,
                 const std::vector<int> &Pos_local_x,
                 const std::vector<int> &Pos_local_y,
                 const std::vector<double> &Number_local,
                 const std::vector<double> &Uold,
                 std::vector<double> &U,
                 const int m);

  int size = F_local.size();
  std::vector<double> x(size);

  std::vector<double> r = F_local;

  double rho = stdvector_dot_mult(r, r); // inner product of two vectors  -> across all processors
  rho = GlobalSum(my_rank, comm_sz, rho);

  // if (my_rank == 0)
  //{
  //   std::cout << "rho=" << rho << "\n";
  // }

  std::vector<double> p = r;

  double residual = stdGetPower(r); // norm of residual -> across all processors
  residual = sqrt(GlobalSum(my_rank, comm_sz, residual));

  // if (my_rank == 0)
  //{
  //   std::cout << "residual=" << residual << "\n";
  // }

  std::vector<double> q(size);
  std::vector<double> temp(size);

  double alpha, rho_old, beta, residual_old, residual_ori, local_dot_prod, global_dot_prod;

  residual_ori = residual;
  for (int iter = 1; iter <= 50000; iter++)
  {
    // q = A*p
    MatMult2D(my_rank, comm_sz, Pos_local_x, Pos_local_y, Number_local, p, q, m); // matrix-vector product -> across all processors

    local_dot_prod = stdvector_dot_mult(p, q); // inner product of two vectors  -> across all processors
    global_dot_prod = GlobalSum(my_rank, comm_sz,
                                local_dot_prod);
    alpha = rho / global_dot_prod;
    // if (my_rank == 0)
    //{
    //   std::cout << "alpha=" << alpha << "\n";
    // }
    temp = scale_stdvector(p, alpha);
    x = stdvector_add(x, temp);
    temp = scale_stdvector(q, -alpha);
    r = stdvector_add(r, temp);
    rho_old = rho;

    rho = stdvector_dot_mult(r, r); // inner product of two vectors -> across all processors
    rho = GlobalSum(my_rank, comm_sz, rho);

    beta = rho / rho_old;
    // if (my_rank == 0)
    //{
    //   std::cout << "beta=" << beta << "\n";
    // }

    temp = scale_stdvector(p, beta);
    p = stdvector_add(r, temp);

    residual_old = residual;
    residual = stdGetPower(r); // norm of residual -> across all processors
    residual = sqrt(GlobalSum(my_rank, comm_sz, residual));

    if (residual / residual_ori < 1e-13)
    {
      break;
    }
  }

  return x;
}

void all_together(const int my_rank,
                  const int comm_sz,
                  const std::vector<double> &pts_position, std::vector<double> &pts_position_all)
{
  int nProc = comm_sz; // be careful
  int numNodes = pts_position.size();
  std::vector<int> eachProcData(nProc);
  MPI_Allgather(&numNodes, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> disp(nProc, 0);
  for (int i = 1; i < disp.size(); i++)
  {
    disp[i] = disp[i - 1] + eachProcData[i - 1];
  }

  int totalProcData = 0;
  for (int i = 0; i < nProc; i++)
  {
    totalProcData += eachProcData[i];
  }

  if (my_rank == 0)
  {
    pts_position_all.resize(totalProcData);
  }

  MPI_Gatherv(pts_position.data(), pts_position.size(), MPI_DOUBLE, pts_position_all.data(), eachProcData.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // MPI_Allgatherv(pts_position.data(), pts_position.size(), MPI_DOUBLE, pts_position_all.data(), eachProcData.data(), disp.data(), MPI_DOUBLE, MPI_COMM_WORLD);
}
