#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
   int N=53;

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
