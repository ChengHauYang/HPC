#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
   int N;

   // read-in N
   printf("\n Input Last Odd Integer:");
   scanf("%i", &N);


   if (N<0 || N%2==0)
   {
      printf(" Invalid value N = %i.\n",N);
      printf(" N must satisfy: N > 0 and N is Odd Number\n");
      exit(1);
   }

   // Output
   int k=0;
   for (int i=1; i<=N; i=i+2)
   {
      k++;
      printf("%i   ", i);
      if (k%5==0 && k !=1){
        printf("\n");
      }
   }

   if (N%10 != 9){
     printf("\n");
   }

   return 0;
}
