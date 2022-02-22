#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fraction.h"

void fraction_reduce(fraction* sum)
{
   void get_prime_factors(int n,
                          int prime_list[],
                          int* num_primes);

   if ( (sum->numer < 0) && (sum->denom < 0) )
   { sum->numer = abs(sum->numer);
     sum->denom = abs(sum->denom); }

   int prime1[100]; int num_prime_1;
   int msign1 = 1; if (sum->numer < 0) { msign1 = -1; }
   sum->numer = abs(sum->numer);
   get_prime_factors(sum->numer,prime1,&num_prime_1);

   int prime2[100]; int num_prime_2;
   int msign2 = 1; if (sum->denom < 0) { msign2 = -1; }
   sum->denom = abs(sum->denom);
   get_prime_factors(sum->denom,prime2,&num_prime_2);

   int  i = 0; int  j = 0;
   int z1 = prime1[i]; int z2 = prime2[j];

   while(i<num_prime_1 && j<num_prime_2)
   {
      if (z1==z2)
      {
         sum->numer = sum->numer/z1;
         sum->denom = sum->denom/z2;

         i  = i+1;
         j  = j+1;
         z1 = prime1[i];
         z2 = prime2[j];
      }
      else
      {
         if (z1>z2)
         {
            j = j+1;
            z2 = prime2[j];
         }
         else
         {
            i = i+1;
            z1 = prime1[i];
         }
      }
   }

   sum->numer = sum->numer*msign1;
   sum->denom = sum->denom*msign2;
}
