#include <TargetConditionals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int main(int argc, char* argv[])
{
   void usage(const char* prog_name);
   //double CompTrap(const double a,const double b,const int N);
   double CompSimps(const double a,const double b,const int N);

/*
   if (argc != 2) { usage(argv[0]); }
   const int N = strtol(argv[1], NULL, 10);
   if (N<1) { usage(argv[0]); }
   assert(N%2==0 && "in Simpson’s Rule you must make N an even number");
*/

   const double a = -1.0; const double b = 1.0;
   //double T = CompTrap(a,b,N);
   double err[8];
   int N[8]={10, 20, 40, 80, 160, 320, 640, 1280};

   const double Iex = 0.62661737464261433833;
   double T;
   for (int i=0;i<8;i++){
     T = CompSimps(a,b,N[i]);
     err[i] = fabs(T-Iex);
//     printf("\n N = %i, T = %23.15e, err = %12.5e\n\n",N[i],T,err);
   }


   printf("====================================================================== \n");
   printf("N\t\t Error\t\t  Ratio of two consecutive errors\n"); //titles
   printf("======================================================================\n");
   for (int i=0;i<8;i++){
     T = CompSimps(a,b,N[i]);
     err[i] = fabs(T-Iex);
     if (i!=7){
       printf("%i\t\t %12.5e\t\t %12.5e\n", N[i], err[i], err[i]/err[i+1]);
     } else{
       printf("%i\t\t %12.5e\n", N[i], err[i]);
     }
   }
   printf("-------------------------------------------------------------------- \n");
   printf("the error divided by %3.2e every time we double N \n",err[6]/err[7]);

   return 0;
}

void usage(const char *prog_name)
{
   fprintf(stderr, "usage: %s <num_intervals>\n", prog_name);
   fprintf(stderr, "   num_intervals should be positive\n");
   exit(1);
}

double CompTrap(const double a, const double b, const int N)
{
   double func(const double x);
   double h = (b-a)/((double)N);
   double T = 0.5*(func(a)+func(b));
   double x = a;

   for (int i=1; i<N; i++)
   {
      x += h;
      T += func(x);
   }

   return h*T;
}

double CompSimps(const double a, const double b, const int N)
{
   double func(const double x);
   double h = (b-a)/((double)N);
   double T = func(a)+func(b)+4*func(a+(N-1)*h);

   for (int i=1; i<N/2; i++)
   {
     T += 2*func(a+2*i*h)+4*func(a-h+2*i*h);
   }

   return h*T/3;
}

double func(const double x)
{
   return exp(-8*x*x);
}
