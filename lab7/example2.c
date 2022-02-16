#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
   float Ntemp;

   // read-in N
   printf("\n Input Last Odd Integer:");
   scanf("%16e", &Ntemp);

   if (fabs(Ntemp-(int)Ntemp)>1e-8)
   {
      printf(" Invalid value N = %16e.\n",Ntemp);
      printf(" N must satisfy: N is Integer, N > 0,  and N is Odd Number\n");
      exit(1);
   }

   int N=(int)Ntemp;

   if (N<0 || N%2==0)
   {
      printf(" Invalid value N = %i.\n",N);
      printf(" N must satisfy: N is Integer, N > 0,  and N is Odd Number\n");
      exit(1);
   }


   // Output
   int k=0;
   for (int i=1; i<=N; i=i+2)
   {
      k++;
      printf("%5i   ", i);
      if (k%5==0 && k !=1){
        printf("\n");
      }
   }

   if (N%10 != 9){
     printf("\n");
   }

   return 0;
}
