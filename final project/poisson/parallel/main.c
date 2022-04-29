#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "math.h"
#include <mpi.h>
#include <vector>
#include <iostream>

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
    // printf("master\n");
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

void print_myrank(const int my_rank)
{
  printf("my_rank = %d\n", my_rank);
  MPI_Barrier(MPI_COMM_WORLD);
}

void print_master(const int my_rank)
{
  if (my_rank == 0)
  {
    printf("master\n");
  }
}

void print_int_vector(std::vector<int> const &input)
{
  for (int i = 0; i < input.size(); i++)
  {
    std::cout << input.at(i) << ' ';
  }
  std::cout << "\n";
  std::cout << "==================================\n";
}

void print_double_vector(const int my_rank, std::vector<double> const &input)
{
  //std::cout<<"my_rank="<<my_rank<<"\n";
  for (int i = 0; i < input.size(); i++)
  {
    std::cout << input.at(i) << ' ';
  }
  std::cout << "\n";
  std::cout << "==================================\n";
}



void print_double_vector_seq(const int comm_sz,const int my_rank, std::vector<double> const &input)
{
  for (int i=0;i<comm_sz;i++){
    if (my_rank==i){
    std::cout<<"my_rank="<<i<<"\n";
    for (int i = 0; i < input.size(); i++)
    {
      std::cout << input.at(i) << ' ';
    }
    std::cout << "\n";
    std::cout << "==================================\n";    
    }
  }
    MPI_Barrier(MPI_COMM_WORLD);

}

int main(int argc, char *argv[])
{

  // Setup MPI code
  int comm_sz, my_rank;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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


  /// Boundary
  if (my_rank == 0)
  {
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
  MPI_Type_contiguous(3, MPI_INT, &VectortypeInt);
  MPI_Type_commit(&VectortypeInt);

  for (int e = 1; e <= Local_Nele; e++)
  { // Loop over element

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
      // printf("sending\n");
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
        // printf("receiving\n");
        MPI_Recv(&Ke_local, 1, Matrixtype, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&Fe_local, 1, Vectortype, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&num_global, 1, VectortypeInt, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        count_num = -1;
        for (int i = 1; i <= 3; i++)
        {
          for (int j = 1; j <= 3; j++)
          {
            count_num++;
            // printf("Ke_local=%23.15f \n",Ke_local[count_num]);
            // printf("num_global=%d \n", num_global[i - 1]);

            mget(K, num_global[i - 1], num_global[j - 1]) = mget(K, num_global[i - 1], num_global[j - 1]) + Ke_local[count_num];
          }
        }

        for (int i = 1; i <= 3; i++)
        {
          vget(F, num_global[i - 1]) = vget(F, num_global[i - 1]) + Fe_local[i - 1];
        }
      }
    }
    //    MPI_Barrier(MPI_COMM_WORLD);
  }

  // print_matrix(&K);
  // print_vector(&F);

  /// solve Linear System
  // vector u = solve(&K,&F);

  int DivideNodesNum_processor = 0;
  int Sum_DivideNodesNum_processor = 0;
  int sparse_size;
  std::vector<int> Pos_local_x = {};
  std::vector<int> Pos_local_y = {};
  std::vector<double> Number_local = {};
  std::vector<double> F_local = {};

  std::vector<int> DivideNodesNum_processor_proc(comm_sz - 1);
  std::vector<int> sparse_size_proc(comm_sz - 1);
  std::vector<std::vector<int>> Pos_local_x_proc(comm_sz - 1);
  std::vector<std::vector<int>> Pos_local_y_proc(comm_sz - 1);
  std::vector<std::vector<double>> Number_local_proc(comm_sz - 1);
  std::vector<std::vector<double>> F_local_proc(comm_sz - 1);

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

    int DivideNodesNum[comm_sz];
    int commun_num = floor(m / comm_sz);

    // how many nodes
    for (int i = 0; i < comm_sz - 1; i++)
    {
      DivideNodesNum[i] = commun_num;
    }
    DivideNodesNum[comm_sz - 1] = m - commun_num * (comm_sz - 1);

    DivideNodesNum_processor = DivideNodesNum[0] * m;
    //printf("DivideNodesNum_processor=%d \n", DivideNodesNum_processor);
    for (int b = 1; b <= DivideNodesNum_processor; b++)
    {
      for (int a = 1; a <= Global_Nnodes; a++)
      {
        if (mget(K, Sum_DivideNodesNum_processor + b, a) != 0)  // meg(K,i,j)*veg(b,j)
        {
          Pos_local_x.push_back(a);
//          Pos_local_y.push_back(Sum_DivideNodesNum_processor + b);
          Pos_local_y.push_back(b);
          Number_local.push_back(mget(K, Sum_DivideNodesNum_processor + b, a));
        }
      }
    }
    //printf("finish K_local\n");

    F_local.resize(DivideNodesNum_processor);
    for (int j = 1; j <= DivideNodesNum_processor; j++)
    {
      F_local[j - 1] = vget(F, Sum_DivideNodesNum_processor + j);
    }
    //printf("finish F_local\n");

    Sum_DivideNodesNum_processor += DivideNodesNum_processor;

    int Sum_DivideNodesNum_processor_before1=0;
    for (int i = 1; i < comm_sz; i++)
    {
      if (i>1){
        Sum_DivideNodesNum_processor_before1 = Sum_DivideNodesNum_processor - DivideNodesNum_processor_proc[i - 2];
      }

      DivideNodesNum_processor_proc[i - 1] = DivideNodesNum[i] * m;
      //printf("DivideNodesNum_processor=%d \n", DivideNodesNum_processor_proc[i - 1]);

      Pos_local_x_proc[i - 1] = {};
      Pos_local_y_proc[i - 1] = {};
      Number_local_proc[i - 1] = {};

      F_local_proc[i - 1] = {};
      F_local_proc[i - 1].resize(DivideNodesNum_processor_proc[i - 1]);

      for (int b = 1; b <= DivideNodesNum_processor_proc[i - 1]; b++)
      {
        for (int a = 1; a <= Global_Nnodes; a++)
        {

          if (mget(K, Sum_DivideNodesNum_processor + b, a) != 0)
          {
            Pos_local_x_proc[i - 1].push_back(a - Sum_DivideNodesNum_processor_before1);
//            Pos_local_y_proc[i - 1].push_back(Sum_DivideNodesNum_processor + b);
            Pos_local_y_proc[i - 1].push_back(b);
            Number_local_proc[i - 1].push_back(mget(K, Sum_DivideNodesNum_processor + b, a));
          }
          // K_local[count_num] = mget(K, a, Sum_DivideNodesNum_processor + b);
        }
      }
      //std::cout << "0. size: " << Number_local_proc[i - 1].size() << '\n';
      //printf("finish K_local\n");

      sparse_size_proc[i - 1] = Number_local_proc[i - 1].size();

      for (int b = 1; b <= DivideNodesNum_processor_proc[i - 1]; b++)
      {
        F_local_proc[i - 1][b - 1] = vget(F, Sum_DivideNodesNum_processor + b);
      }
      //printf("finish F_local\n");

      Sum_DivideNodesNum_processor += DivideNodesNum_processor_proc[i - 1];

      MPI_Send(&DivideNodesNum_processor_proc[i - 1], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&sparse_size_proc[i - 1], 1, MPI_INT, i, 0, MPI_COMM_WORLD);

      MPI_Send(&F_local_proc[i - 1][0], DivideNodesNum_processor_proc[i - 1], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      //printf("sending vector\n");
      MPI_Send(&Pos_local_x_proc[i - 1][0], sparse_size_proc[i - 1], MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&Pos_local_y_proc[i - 1][0], sparse_size_proc[i - 1], MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&Number_local_proc[i - 1][0], sparse_size_proc[i - 1], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

      //printf("sending matrix\n");
    }
  }
  else
  {

    MPI_Recv(&DivideNodesNum_processor, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&sparse_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //std::cout << "receiving number\n";
    F_local.resize(DivideNodesNum_processor);
    Pos_local_x.resize(sparse_size);
    Pos_local_y.resize(sparse_size);
    Number_local.resize(sparse_size);

    MPI_Recv(&F_local[0], DivideNodesNum_processor, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //std::cout << "receiving F\n";

    MPI_Recv(&Pos_local_x[0], sparse_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //std::cout << "receiving pos_x\n";

    MPI_Recv(&Pos_local_y[0], sparse_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //std::cout << "receiving pos_y\n";

    MPI_Recv(&Number_local[0], sparse_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //std::cout << "receiving Number_local\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<double> solveCGMPI(const int my_rank, const int comm_sz,
                                 const std::vector<double> &F_local,
                                 const std::vector<int> &Pos_local_x,
                                 const std::vector<int> &Pos_local_y,
                                 const std::vector<double> &Number_local);
  std::vector<double> u = solveCGMPI(my_rank, comm_sz, F_local, Pos_local_x, Pos_local_y, Number_local);

  void all_together(const int my_rank,
                    const int comm_sz,
                    const std::vector<double> &pts_position, std::vector<double> &pts_position_all);

  std::vector<double> u_all;

  //print_double_vector(my_rank,u);

  all_together(my_rank, comm_sz, u, u_all);

  if (my_rank == 0)
  {
    //print_double_vector(my_rank,u_all);
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
      fprintf(outfile, "%25.20e\n", u_all[i - 1]);
    }

    // close file
    fclose(outfile);

    // Call python script to plot
    system("python3.8 phase_plot.py");
  }

  /*

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

    MPI_Barrier(MPI_COMM_WORLD);
    */
}
