#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

int main(int argc, char* argv[])
{
   void usage(const char* prog_name);
   double CompTrap(const double a,const double b,
		   const int N,const int thread_count);
   double CompSimps(const double a,const double b,
		   const int N,const int thread_count);

   if (argc != 3) { usage(argv[0]); }
   const int thread_count = strtol(argv[1], NULL, 10);
   const int N = strtol(argv[2], NULL, 10);
   if (thread_count<1 || N<1 || N % thread_count != 0)
   { usage(argv[0]); }

   const double a = -1.0; const double b = 1.0;
   const double time1 = omp_get_wtime();
   //double T = CompTrap(a,b,N,thread_count);
   double T = CompSimps(a,b,N,thread_count);
   const double time2 = omp_get_wtime();
   const double Iex = 0.62661737464261433833;
   double err = fabs(T-Iex);

   printf("\n N = %i, T = %23.15e, err = %12.5e,\n",N,T,err);
   printf(" time = %12.5e\n\n",time2-time1);

   return 0;
}

void usage(const char *prog_name)
{
   fprintf(stderr,"usage: %s <num_threads> <num_intervals>\n",
	   prog_name);
   fprintf(stderr,"   num_threads should be positive\n");
   fprintf(stderr,"   num_intervals should be positive\n");
   fprintf(stderr,"   mod(num_intervals,num_threads) != 0\n");
   exit(1);
}

double CompTrap(const double a, const double b,
		const int N, const int thread_count)
{
   double func(const double x);
   double h = (b-a)/((double)N);

   double T = 0.5*(func(a)+func(b));
#  pragma omp parallel for num_threads(thread_count)	\
   reduction(+: T)
   for (int i=1; i<N; i++)
   {
      T += func(a+i*h);
   }
   printf("%d",N);

   return h*T;
}

double CompSimps(const double a, const double b,
		const int N, const int thread_count)
{
   double func(const double x);
   double h = (b-a)/((double)N);

   double T = func(a)+func(b)+4*func(a+(N-1)*h);
#  pragma omp parallel for num_threads(thread_count)	\
   reduction(+: T)
   for (int i=1; i<N/2; i++)
   {
     T += 2*func(a+2*i*h)+4*func(a-h+2*i*h);
   }
   //printf("%d",N);

   return h*T/3;
}

double func(const double x)
{
   return exp(-8*x*x);
}
