#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#include "vector.h"

int main(int argc, char* argv[])
{
  // Functions needed
  void WriteSoln(const int my_rank,
                 const vector* x,
                 const vector* U,
                 const vector* V,
                 const double time,
                 const char* filename);
  double func(const double x);
  void send_boundary_data(const int my_rank,
                          const int comm_sz,
                          const double U1,
                          const double UN);
  void receive_boundary_data(const int my_rank,
                             const int comm_sz,
                             const double alpha,
                             const double beta,
                             double* alpha_local,
                             double* beta_local);

  // Setup MPI code
  int comm_sz, my_rank;
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  FILE* file = fopen("num_processors.data","w");
  fprintf(file,"%i\n",comm_sz);
  fclose(file);

  // Grab the global N parameter
  // and set the local N parameter
  int N;
  void get_input(int argc, char* argv[],
                 const int my_rank,
                 const int comm_sz,
                 int* N);
  get_input(argc,argv,my_rank,comm_sz,&N);
  const int N_local = N/comm_sz;
  // must have at least 1 interior node per processor
  assert(N_local>=1);
  const double h = 1.0/((double)N);

  // Dimension vectors
  vector Unew = new_vector(N_local);
  vector U    = new_vector(N_local);
  vector Vnew = new_vector(N_local);  
  vector V    = new_vector(N_local);
  vector x    = new_vector(N_local);
  double time = 0.0;

  // Initial condition
  for (int i=1; i<=N_local; i++)
    {
      double xtmp = ((double)(i + my_rank*N_local)-0.5)*h;
      vget(x,i) = xtmp;
      vget(U,i) = func(xtmp);
      vget(V,i) = func(xtmp);
    }

  // Print initial guess to a file
  char proc_number[10];
  char frame_number[10];
  char filename[80];
  sprintf(proc_number,"%i",my_rank);
  sprintf(frame_number,"%i",0);
  strcpy(filename, "q");
  strcat(filename, frame_number);
  strcat(filename, "_");
  strcat(filename, proc_number);
  strcat(filename, ".data");
  WriteSoln(my_rank,&x,&U,&V,time,filename);

  // Get time
  double time_start;
  if (my_rank==0)
    { time_start = MPI_Wtime(); }

  // -----------------------------------
  // Main time stepping loop
  // -----------------------------------
  double dt = 0.95*h;
  const double Tfinal = 1.0;
  const int NumFrames = 10;
  const double T_per_frame = Tfinal/((double)NumFrames);
  const int NumSteps_per_frame = (int)(T_per_frame/dt+0.5);
  dt = T_per_frame/((double)NumSteps_per_frame); // reset dt
  double nu = dt/h; double nu2 = pow(nu,2); // CFL number
  const int NumSteps = NumSteps_per_frame*NumFrames;

  for (int nframe=1; nframe<=NumFrames; nframe++)
    {
      for (int n=1; n<=NumSteps_per_frame; n++)
      {
      // send and receive boundary data
      send_boundary_data(my_rank,comm_sz,
                         vget(U,1),vget(U,N_local));
      double ULBC,URBC;
      receive_boundary_data(my_rank,comm_sz,
                            vget(U,N_local),vget(U,1),
                            &ULBC,&URBC);

      send_boundary_data(my_rank,comm_sz,
                         vget(V,1),vget(V,N_local));
      double VLBC,VRBC;
      receive_boundary_data(my_rank,comm_sz,
                            vget(U,N_local),vget(V,1),
                            &VLBC,&VRBC);

      // time-stepping method for U
      time += dt;
      vget(Unew,1) = vget(U,1) - 0.5*nu*(vget(U,2)-ULBC)
                   + 0.5*nu2*(vget(U,2)-2.0*vget(U,1)+ULBC);
      for (int i=2; i<N_local; i++)
        {
          vget(Unew,i) = vget(U,i) - 0.5*nu*(vget(U,i+1)-vget(U,i-1))
                       + 0.5*nu2*(vget(U,i+1)-2.0*vget(U,i)+vget(U,i-1));
        }
      vget(Unew,N_local) = vget(U,N_local) - 0.5*nu*(URBC-vget(U,N_local-1))
                         + 0.5*nu2*(URBC-2.0*vget(U,N_local)+vget(U,N_local-1));

      // time-stepping method for V
      vget(Vnew,1) = vget(V,1) + 0.5*nu*(vget(V,2)-VLBC)
                   + 0.5*nu2*(vget(V,2)-2.0*vget(V,1)+VLBC);
      for (int i=2; i<N_local; i++)
        {
          vget(Vnew,i) = vget(V,i) + 0.5*nu*(vget(V,i+1)-vget(V,i-1))
                       + 0.5*nu2*(vget(V,i+1)-2.0*vget(V,i)+vget(V,i-1));
        }
      vget(Vnew,N_local) = vget(V,N_local) - 0.5*nu*(VRBC-vget(V,N_local-1))
                         + 0.5*nu2*(VRBC-2.0*vget(V,N_local)+vget(V,N_local-1));

      // reset for next step
      copy_vector(&Unew,&U);
      copy_vector(&Vnew,&V);
      }

      // Print final guess to a file
      sprintf(frame_number,"%i",nframe);
      strcpy(filename, "q");
      strcat(filename, frame_number);
      strcat(filename, "_");
      strcat(filename, proc_number);
      strcat(filename, ".data");
      WriteSoln(my_rank,&x,&U,&V,time,filename);
    }

  // Print elapsed time to screen
  if (my_rank==0)
    {
      double time_end = MPI_Wtime();
      double time_elapsed = time_end - time_start;
      printf(" NP = %2i, N = %i, dt = %11.5e, nu = %11.5e, NumSteps = %i, Elapsed time = %20.13e\n",
             comm_sz,N,dt,nu,NumSteps,time_elapsed);
    }

  // Free up memory
  delete_vector(&Unew);
  delete_vector(&U);
  delete_vector(&V);
  delete_vector(&x);

  // End program
  MPI_Finalize();
  return 0;
}

void get_input(int argc, char* argv[],
               const int my_rank,
               const int comm_sz,
               int* N)
{
  void usage(const char *prog_name);
  if (my_rank==0)
    {
      if (argc != 2) { usage(argv[0]); }

      *N = strtol(argv[1], NULL, 10);
      if (*N<=0) { usage(argv[0]); }
      if (*N%comm_sz!=0) { usage(argv[0]); }

      for (int i=1; i<comm_sz; i++)
        {
          MPI_Send(N,1,MPI_INT,   i,0,MPI_COMM_WORLD);
        }
    }
  else
    {
      MPI_Recv(N,1,MPI_INT,0,0,MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
}

void usage(const char *prog_name)
{
  fprintf(stderr,"usage: %s <N>\n",prog_name);
  fprintf(stderr,"   N should be positive\n");
  fprintf(stderr,"   N should be exactly divisible "
          "by the number of processors\n");
  exit(1);
}

void WriteSoln(const int my_rank,
               const vector* x,
               const vector* U,
               const vector* V,
               const double time,
               const char* filename)
{
  const int N = U->size;
  FILE* file = fopen(filename,"w");

  fprintf(file,"%24.15e\n",time);
  for (int i=1; i<=N; i++)
    {
      const double xtmp = vgetp(x,i);
      const double Utmp = vgetp(U,i);
      const double Vtmp = vgetp(V,i);
      fprintf(file,"%24.15e %24.15e %24.15e\n",xtmp,Utmp,Vtmp);
    }
  fclose(file);
}

void send_boundary_data(const int my_rank,
                        const int comm_sz,
                        const double U1,
                        const double UN)
{
  const int last_rank = comm_sz-1;
  if (comm_sz==1)
    { return; }
  else
    {
      if (my_rank==0)
        {
          MPI_Send(&UN,1,MPI_DOUBLE,1,0,
                   MPI_COMM_WORLD);
          MPI_Send(&U1,1,MPI_DOUBLE,last_rank,0,
                   MPI_COMM_WORLD);
        }
      else if (my_rank==last_rank)
        {
          MPI_Send(&UN,1,MPI_DOUBLE,0,0,
                   MPI_COMM_WORLD);
          MPI_Send(&U1,1,MPI_DOUBLE,last_rank-1,0,
                   MPI_COMM_WORLD);
        }
      else
        {
          MPI_Send(&U1,1,MPI_DOUBLE,my_rank-1,0,
                   MPI_COMM_WORLD);
          MPI_Send(&UN,1,MPI_DOUBLE,my_rank+1,0,
                   MPI_COMM_WORLD);
        }
    }
}

void receive_boundary_data(const int my_rank,
                           const int comm_sz,
                           const double Uleft,
                           const double Uright,
                           double* Uleft_local,
                           double* Uright_local)
{
  const int last_rank = comm_sz-1;
  if (comm_sz==1)
    {
      *Uleft_local  = Uleft;
      *Uright_local = Uright;
    }
  else
    {
      // THEN RECEIVE
      if (my_rank==0)
        {
          MPI_Recv(Uleft_local,1,MPI_DOUBLE,last_rank,0,
                   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          MPI_Recv(Uright_local,1,MPI_DOUBLE,1,0,
                   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      else if (my_rank==last_rank)
        {
          MPI_Recv(Uleft_local,1,MPI_DOUBLE,last_rank-1,0,
                   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          MPI_Recv(Uright_local,1,MPI_DOUBLE,0,0,
                   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      else
        {
          MPI_Recv(Uleft_local,1,MPI_DOUBLE,my_rank-1,0,
                   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          MPI_Recv(Uright_local,1,MPI_DOUBLE,my_rank+1,0,
                   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
}

double func(const double x)
{
  return exp(-400.0*pow(x-0.5,2));
}
