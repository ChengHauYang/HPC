#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "vector.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char *argv[])
{
  // Setup MPI code
  int comm_sz, my_rank;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // Grab the global N parameter
  // and set the local N parameter

  int num_test = 6;
  int Nall[num_test];
  double times[num_test];
  double pi_cal[num_test];
  for (int i = 0; i < num_test; i++)
  {
    Nall[i] = pow(10, i + 1);
  }

  for (int i = 0; i < num_test; i++)
  {
    int N = Nall[i];
    const int N_local = N / comm_sz;

    // Get time
    double time_start;
    if (my_rank == 0)
    {
      time_start = MPI_Wtime();
    }

    /*
      // Create partial vector on current processor
      vector v_local = new_vector(N_local);
      for (int i = 1; i <= N_local; i++)
      {
        vget(v_local, i) = sqrt((double)(i + my_rank * N_local));
      }

      // Compute 2-norm squared
      double norm_squared;
      double norm_squared_local = pow(vget(v_local, 1), 2);
      for (int i = 2; i <= N_local; i++)
      {
        norm_squared_local += pow(vget(v_local, i), 2);
      }

      // Add local results to the global result on Processor 0
      if (my_rank != 0)
      {
        MPI_Send(&norm_squared_local, 1, MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
      }
      else
      {
        norm_squared = norm_squared_local;
        for (int i = 1; i < comm_sz; i++)
        {
          MPI_Recv(&norm_squared_local, 1, MPI_DOUBLE, i, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          norm_squared += norm_squared_local;
        }
      }
    */

    double x, y, r, pi;
    int in_circle;
    int in_circle_local;
    in_circle_local = 0;
    // printf("N_local=%i\n",N_local);
    for (int k = 0; k <= N_local; k++)
    {
      x = ((double)rand() / (RAND_MAX));
      y = ((double)rand() / (RAND_MAX));
      r = sqrt(x * x + y * y);
      if (r <= 1)
      {
        in_circle_local++;
      }
    }

    // Add local results to the global result on Processor 0
    if (my_rank != 0)
    {
      MPI_Send(&in_circle_local, 1, MPI_INT, 0, 0,
               MPI_COMM_WORLD);
    }
    else
    {
      in_circle = in_circle_local; // master one
      for (int i = 1; i < comm_sz; i++)
      {
        MPI_Recv(&in_circle_local, 1, MPI_INT, i, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // printf("in_circle_local=%i\n",in_circle_local);
        in_circle += in_circle_local;
      }
      pi = (double)4 * in_circle / N;
      pi_cal[i] = pi;
    }

    // Print answer to screen
    if (my_rank == 0)
    {
      double time_end = MPI_Wtime();
      double time_elapsed = time_end - time_start;
      printf(" NP = %2i, N = %i, relative error = %20.13e\n",
             comm_sz, N, fabs((pi - M_PI) / M_PI));
      printf("     Elapsed time = %20.13e\n",
             time_elapsed);
      times[i] = time_elapsed;
    }
  }
  if (my_rank == 0)
  {

    /// Print solution to file
    char filename[] = "mc_pi_table.txt";

    // open file
    FILE *outfile = fopen(filename, "w");

    // output data
    fprintf(outfile, "N\t\t,cpu time\t\t,relative error \n");
    for (int i = 0; i < num_test; i++)
    {
      fprintf(outfile, "%d\t\t", Nall[i]);
      fprintf(outfile, ",%25.20e\t\t", times[i]);
      fprintf(outfile, ",%25.20e\n", fabs((pi_cal[i] - M_PI) / M_PI));
    }

    // close file
    fclose(outfile);
  }
  // End program
  MPI_Finalize();
  return 0;
}
