#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#ifdef _OPENMP
#  include <omp.h>
#endif

int main()
{

   const int thread_count_max = 16;
   int num_test = 6;
   int N[num_test];
   double times[num_test];
   double pi_cal[num_test];
   for (int i=0;i<num_test;i++){
     N[i]=pow(10,i+1  );
   }

   double x,y,r,pi;
   int in_circle;
   double time1,time2,clock_time;
   for (int i=0;i<num_test;i++){
     time1 = omp_get_wtime();
     in_circle=0;
# pragma omp parallel for num_threads(thread_count_max) \
     reduction(+: in_circle)
     for(int k=0;k<=N[i];k++){
       x =  ((double) rand() / (RAND_MAX));
       y =  ((double) rand() / (RAND_MAX));
       r = sqrt(x*x+y*y);
       if (r<=1){
         in_circle++;
       }
     }

     pi = (double)4*in_circle/N[i];
     pi_cal[i]=pi;
     time2 = omp_get_wtime();
     clock_time = time2-time1;
     times[i]=clock_time;

   }



   printf("====================================================================== \n");
   printf("N\t\t Cpu Time \t\t  Relative Error\n"); //titles
   printf("======================================================================\n");
   for (int i=0;i<num_test;i++){
     printf("%i\t\t %13.6e\t\t %13.6e\n", N[i], times[i], fabs((pi_cal[i]-M_PI)/M_PI));
   }
   printf("-------------------------------------------------------------------- \n");

   /// Print solution to file
   char filename[] = "mc_pi_table.txt";

   // open file
   FILE* outfile = fopen(filename,"w");

   // output data
   fprintf(outfile,"N\t\t,cpu time\t\t,relative error \n");
   for (int i=0;i<num_test;i++){
     fprintf(outfile,"%d\t\t",N[i]);
     fprintf(outfile,",%25.20e\t\t",times[i]);
     fprintf(outfile,",%25.20e\n",fabs((pi_cal[i]-M_PI)/M_PI));
   }

   // close file
   fclose(outfile);

   return 0;
}
