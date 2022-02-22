#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void get_prime_factors(int n,
                       int prime_list[],
                       int* num_primes)
{
   *num_primes = 0;

   while (n%2==0)
   {
      prime_list[*num_primes] = 2;
      *num_primes = *num_primes + 1;
      n = n/2;
   }

   for (int i=3; i<=sqrt(n); i=i+2)
   {
      while (n%i==0)
      {
         prime_list[*num_primes] = i;
         *num_primes = *num_primes + 1;
         n = n/i;
      }
   }

   if (n>2)
   {
      prime_list[*num_primes] = n;
      *num_primes = *num_primes + 1;
   }
}
