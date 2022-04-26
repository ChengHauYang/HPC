#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "math.h"
#include <mpi.h>

void usage(const char *prog_name)
{
  fprintf(stderr, "usage: %s <N>\n", prog_name);
  fprintf(stderr, "   N should be positive\n");
  fprintf(stderr, "   N-1 should be exactly divisible "
                  "by the number of processors\n");
  exit(1);
}

void get_input(int argc, char *argv[],
               const int my_rank,
               const int comm_sz,
               int *N)
{
  void usage(const char *prog_name);
  if (my_rank == 0)
  {
    if (argc != 2)
    {
      usage(argv[0]);
    }

    *N = strtol(argv[1], NULL, 10);
    if (*N <= 0)
    {
      usage(argv[0]);
    }
    if ((*N - 1) % comm_sz != 0)
    {
      usage(argv[0]);
    }

    for (int i = 1; i < comm_sz; i++)
    {
      MPI_Send(N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
  else
  {
    MPI_Recv(N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }
}

int main(int argc, char *argv[])
{

  // Setup MPI code
  int comm_sz, my_rank;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /*
    int m=0;
    printf("    Input m: ");
    scanf("%i",&m);
  */

  // pick m as 17 is better !
  // m-1 should be exactly divisible by the number of processors

  int m;
  void get_input(int argc, char *argv[],
                 const int my_rank,
                 const int comm_sz,
                 int *m);
  get_input(argc, argv, my_rank, comm_sz, &m);

  matrix Global_Coords = new_matrix(m * m, 2); // x,y
  matrix Global_NodesNum = new_matrix(2 * (m - 1) * (m - 1), 3);

  int Global_Nnodes = m * m;
  int Global_Nele = 2 * (m - 1) * (m - 1);
  int Local_Nele = 2 * (m - 1) * (m - 1) / comm_sz;

  matrix Local_NodesNum = new_matrix(2 * (m - 1) * (m - 1) / comm_sz, 3);
  vector boundary = new_vector(4 * m - 4);

  // printf("comm_sz=%d\n",comm_sz);
  // printf("my_rank=%d\n",my_rank);

  /// Coords from 1 ~ m*m
  /*
  int k,i,j;
  for (int num=my_rank*m*m/comm_sz;num<(my_rank+1)*m*m/comm_sz;num++){
      k=num-my_rank*m*m/comm_sz+1;
      j = floor(num/m)+1;
      i = num%m+1;
      mget(Coords, k, 1) = (double)(j-1)/(m-1);
      mget(Coords, k, 2) = (double)(i-1)/(m-1);
  }
  */

  int k = 0;
  for (int i = 1; i <= m; i++)
  {
    for (int j = 1; j <= m; j++)
    {
      k++;
      mget(Global_Coords, k, 1) = (double)(j - 1) / (m - 1);
      mget(Global_Coords, k, 2) = (double)(i - 1) / (m - 1);
    }
  }

  /// Nodes number for every triangle
  int true_num = 1;
  for (int i = 1; i <= m * m - (m + 2) + 1; i = i + 1)
  {
    // printf("%d\n",i);
    if (i % m != 0)
    {
      mget(Global_NodesNum, true_num, 1) = 1 + i - 1;
      mget(Global_NodesNum, true_num, 2) = 2 + i - 1;
      mget(Global_NodesNum, true_num, 3) = m + 2 + i - 1;
      mget(Global_NodesNum, true_num + 1, 1) = 1 + i - 1;
      mget(Global_NodesNum, true_num + 1, 2) = m + 2 + i - 1;
      mget(Global_NodesNum, true_num + 1, 3) = m + 1 + i - 1;
      true_num = true_num + 2;
    }
  }

  // print_matrix(&Global_NodesNum);

  for (int i = 1; i <= Local_Nele; i = i + 1)
  {
    for (int j = 1; j <= 3; j = j + 1)
    {
      mget(Local_NodesNum, i, j) = mget(Global_NodesNum, i + Local_Nele * my_rank, j);
    }
  }

  // if (my_rank==0){
  //   print_matrix(&Local_NodesNum);
  // }

  // printf("my rank=%d", my_rank);
  // print_matrix(&Local_NodesNum);

  /// Boundary
  k = 0;

  for (int i = 1; i <= m; i++)
  { // y-
    k++;
    vget(boundary, k) = i;
  }

  for (int i = m + 1; i <= m * m; i = i + m)
  { // x-
    // printf("%d\n",i);
    k++;
    vget(boundary, k) = i;
  }

  for (int i = 2 * m; i <= m * m; i = i + m)
  { // x+
    k++;
    vget(boundary, k) = i;
  }

  for (int i = m * m - m + 2; i <= m * m - 1; i = i + 1)
  { // y+
    k++;
    vget(boundary, k) = i;
  }
  // print_vector(&boundary);

  /// Stiffness matrix and force vector
  matrix K = new_matrix(Global_Nnodes, Global_Nnodes);
  vector F = new_vector(Global_Nnodes);

  /// Loop over element to do the assembly
  vector tri_nodes_number = new_vector(3);
  matrix Me = new_matrix(3, 3);
  matrix Be = new_matrix(2, 3);
  matrix Ke = new_matrix(3, 3);
  vector Fe = new_vector(3);

  // My own mpi type
  MPI_Datatype Matrixtype;
  MPI_Type_contiguous(9, MPI_DOUBLE, &Matrixtype);
  MPI_Type_commit(&Matrixtype);

  MPI_Datatype Vectortype;
  MPI_Type_contiguous(3, MPI_DOUBLE, &Vectortype);
  MPI_Type_commit(&Vectortype);

  MPI_Datatype VectortypeInt;
  MPI_Type_contiguous(3, MPI_DOUBLE, &VectortypeInt);
  MPI_Type_commit(&VectortypeInt);

  if (my_rank == 0)
  {
    printf("master\n");
  }
  for (int e = 1; e <= Local_Nele; e++)
  { // Loop over element

    if (my_rank == 0)
    {
      printf("master\n");
    }
    int num_global[3];
    for (int i = 1; i <= 3; i++)
    {
      num_global[i - 1] = mget(Local_NodesNum, e, i);
      vget(tri_nodes_number, i) = mget(Local_NodesNum, e, i);
    }
    // print_vector(&tri_nodes_number);
    for (int i = 1; i <= 3; i++)
    {
      int NodesOnTri = vget(tri_nodes_number, i);
      mget(Me, i, 1) = 1;
      mget(Me, i, 2) = mget(Global_Coords, NodesOnTri, 1); // x
      mget(Me, i, 3) = mget(Global_Coords, NodesOnTri, 2); // y
    }
    // print_matrix(&Me);

    double det_33(const matrix *A);
    double Area = det_33(&Me) / 2;
    // printf("Area=%10.9e\n",Area);

    // Be = Element shape function derivative matrices
    mget(Be, 1, 1) = (mget(Global_Coords, (int)vget(tri_nodes_number, 2), 2) - mget(Global_Coords, (int)vget(tri_nodes_number, 3), 2)) / 2 / Area;
    mget(Be, 1, 2) = (mget(Global_Coords, (int)vget(tri_nodes_number, 3), 2) - mget(Global_Coords, (int)vget(tri_nodes_number, 1), 2)) / 2 / Area;
    mget(Be, 1, 3) = (mget(Global_Coords, (int)vget(tri_nodes_number, 1), 2) - mget(Global_Coords, (int)vget(tri_nodes_number, 2), 2)) / 2 / Area;
    mget(Be, 2, 1) = (mget(Global_Coords, (int)vget(tri_nodes_number, 3), 1) - mget(Global_Coords, (int)vget(tri_nodes_number, 2), 1)) / 2 / Area;
    mget(Be, 2, 2) = (mget(Global_Coords, (int)vget(tri_nodes_number, 1), 1) - mget(Global_Coords, (int)vget(tri_nodes_number, 3), 1)) / 2 / Area;
    mget(Be, 2, 3) = (mget(Global_Coords, (int)vget(tri_nodes_number, 2), 1) - mget(Global_Coords, (int)vget(tri_nodes_number, 1), 1)) / 2 / Area;
    // print_matrix(&Be);

    // build local stiffness matrix
    matrix matrix_mult_transpose(const matrix *A, const matrix *B);
    matrix temp_Ke = matrix_mult_transpose(&Be, &Be);
    matrix scale_matrix(const matrix *A, double scaling);
    Ke = scale_matrix(&temp_Ke, Area);
    // print_matrix(&Ke);

    // build local force vector
    // center point
    double xg = (mget(Global_Coords, (int)vget(tri_nodes_number, 1), 1) + mget(Global_Coords, (int)vget(tri_nodes_number, 2), 1) + mget(Global_Coords, (int)vget(tri_nodes_number, 3), 1)) / 3;
    double yg = (mget(Global_Coords, (int)vget(tri_nodes_number, 1), 2) + mget(Global_Coords, (int)vget(tri_nodes_number, 2), 2) + mget(Global_Coords, (int)vget(tri_nodes_number, 3), 2)) / 3;
    double eval_f = 2 * M_PI * M_PI * sin(M_PI * xg) * sin(M_PI * yg);

    for (int i = 1; i <= 3; i++)
    {
      vget(Fe, i) = Area * eval_f / 3;
    }

    // assemby for stiffness
    double Ke_local[9];
    int count_num = -1;
    for (int i = 1; i <= 3; i++)
    {
      for (int j = 1; j <= 3; j++)
      {
        count_num++;
        // printf("i=%d \n",(int)vget(tri_nodes_number,i));
        // printf("j=%d \n",(int)vget(tri_nodes_number,j));
        // mget(K,(int)vget(tri_nodes_number,i),(int)vget(tri_nodes_number,j))
        //=mget(K,(int)vget(tri_nodes_number,i),(int)vget(tri_nodes_number,j))+mget(Ke,i,j);
        Ke_local[count_num] = mget(Ke, i, j);
      }
    }

    // assemby for force
    double Fe_local[3];
    for (int i = 1; i <= 3; i++)
    {
      Fe_local[i - 1] = vget(Fe, i);
    }

    if (my_rank != 0)
    {
      printf("sending\n");
      MPI_Send(&Ke_local, 1, Matrixtype, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&Fe_local, 1, Vectortype, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&num_global, 1, VectortypeInt, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
      // assemby for stiffness
      for (int i = 1; i <= 3; i++)
      {
        for (int j = 1; j <= 3; j++)
        {
          mget(K, (int)vget(tri_nodes_number, i), (int)vget(tri_nodes_number, j)) = mget(K, (int)vget(tri_nodes_number, i), (int)vget(tri_nodes_number, j)) + mget(Ke, i, j);
        }
      }

      // assemby for force
      for (int i = 1; i <= 3; i++)
      {
        vget(F, (int)vget(tri_nodes_number, i)) = vget(F, (int)vget(tri_nodes_number, i)) + vget(Fe, i);
      }

      for (int k = 1; k < comm_sz; k++) // Loop over processors
      {
        printf("receiving\n");
        MPI_Recv(&Ke_local, 1, Matrixtype, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&Fe_local, 1, Vectortype, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&num_global, 1, VectortypeInt, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        count_num = -1;
        for (int i = 1; i <= 3; i++)
        {
          for (int j = 1; j <= 3; j++)
          {
            count_num++;
            mget(K, num_global[i - 1], num_global[j - 1]) = mget(K, num_global[i - 1], num_global[j - 1]) + Ke_local[count_num];
          }
        }

        for (int i = 1; i <= 3; i++)
        {
          vget(F, num_global[i - 1]) = vget(F, num_global[i - 1]) + Fe_local[i - 1];
        }
      }
    }
  }

  // print_matrix(&K);
  // print_vector(&F);

  /// solve Linear System
  // vector u = solve(&K,&F);
  if (my_rank == 0)
  {

    // strong enforce BC => set LHS = 1 and set RHS = 0
    for (int i = 1; i <= 4 * m - 4; i++)
    {
      vget(F, (int)vget(boundary, i)) = 0;
      for (int j = 1; j < m * m; j++)
      {
        mget(K, (int)vget(boundary, i), j) = 0;
        mget(K, j, (int)vget(boundary, i)) = 0;
      }
      mget(K, (int)vget(boundary, i), (int)vget(boundary, i)) = 1;
    }

    print_matrix(&K);

    // TO DO: try to save stiffness matrix as sparse matrix!!!!!!!!!!!!!!!
    vector solveCG(const matrix *A, const vector *b);
    vector u = solveCG(&K, &F);
    // print_vector(&u);

    /// postprosessing

    /// Print solution to file
    char filename[] = "output.data";

    // open file
    FILE *outfile = fopen(filename, "w");

    // output data
    for (int i = 1; i <= m * m; i++)
    {
      fprintf(outfile, "%25.20e  ", mget(Global_Coords, i, 1));
      fprintf(outfile, "%25.20e  ", mget(Global_Coords, i, 2));
      fprintf(outfile, "%25.20e\n", vget(u, i));
    }

    // close file
    fclose(outfile);

    // Call python script to plot
    system("python3.8 phase_plot.py");

    /// Output Tecplot
    char filenametec[] = "output.tec";
    FILE *outfiletec = fopen(filenametec, "w");
    fprintf(outfiletec, "%s \n", "Titile = 'Poisson' ");
    fprintf(outfiletec, "%s \n", "VARIABLES = X, Y, U, Error");
    fprintf(outfiletec, "zone N= %d, E = %d", Global_Nnodes, Global_Nele);
    fprintf(outfiletec, "\n");
    fprintf(outfiletec, "DATAPACKING=POINT, ZONETYPE=FETRIANGLE");
    fprintf(outfiletec, "\n");
    for (int i = 1; i <= Global_Nnodes; i++)
    {
      double Exact = sin(M_PI * mget(Global_Coords, i, 1)) * sin(M_PI * mget(Global_Coords, i, 2));
      fprintf(outfiletec, "%25.20e  %25.20e  %25.20e  %25.20e\n", mget(Global_Coords, i, 1), mget(Global_Coords, i, 2), vget(u, i), fabs(Exact - vget(u, i)));
    }
    for (int i = 1; i <= Global_Nele; i++)
    {
      fprintf(outfiletec, "%d  %d  %d", (int)mget(Global_NodesNum, i, 1), (int)mget(Global_NodesNum, i, 2), (int)mget(Global_NodesNum, i, 3));
      if (i != Global_Nele)
      {
        fprintf(outfiletec, "\n");
      }
    }
  }
}
