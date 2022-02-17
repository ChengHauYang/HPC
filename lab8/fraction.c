#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct fraction fraction;
typedef struct properfraction properfraction; //struct tag

struct fraction
{
   int numer; int denom;
};

struct properfraction
{
   int whole;
   int numer;
   int denom;
};

int main()
{
   fraction a,b,sum,diff,prod,quot;
/*
   a.numer = -1; a.denom = 10;
   b.numer = -3; b.denom = 8;
*/

   properfraction psum,pdiff,pprod,pquot;

   // read-in
   printf("\n Input numerator of a:");
   scanf("%i", &a.numer);
   printf("\n Input denominator of a:");
   scanf("%i", &a.denom);
   printf("\n Input numerator of b:");
   scanf("%i", &b.numer);
   printf("\n Input denominator of b:");
   scanf("%i", &b.denom);

   void fraction_add(const fraction* a,
                     const fraction* b,
                     fraction* sum);
   fraction_add(&a,&b,&sum);

   void fraction_substract(const fraction* a,
                     const fraction* b,
                     fraction* diff);
   fraction_substract(&a,&b,&diff);

   void fraction_mult(const fraction* a,
                     const fraction* b,
                     fraction* prod);
   fraction_mult(&a,&b,&prod);

   void fraction_divide(const fraction* a,
                     const fraction* b,
                     fraction* quot);
   fraction_divide(&a,&b,&quot);

   void make_proper(const fraction* in,
                         properfraction* out); // avoid implicit declaration

   printf("========================================");

   printf("\n %i/%i  +  %i/%i  =  %i/%i\n\n",
          a.numer,a.denom,b.numer,b.denom,
          sum.numer,sum.denom);
   make_proper(&sum,&psum);
   printf("Proper Fraction:");
   if (abs(psum.whole)>0)
     { printf("%i(%i/%i)",psum.whole,psum.numer,psum.denom); }
   else
     { printf("%i/%i",psum.numer,psum.denom); }
   printf("\n");

   printf("\n %i/%i  -  %i/%i  =  %i/%i\n\n",
          a.numer,a.denom,b.numer,b.denom,
          diff.numer,diff.denom);
   make_proper(&diff,&pdiff);
   printf("Proper Fraction:");
   if (abs(pdiff.whole)>0)
    { printf("%i(%i/%i)",pdiff.whole,pdiff.numer,pdiff.denom); }
   else
    { printf("%i/%i",pdiff.numer,pdiff.denom); }
  printf("\n");

   printf("\n %i/%i  *  %i/%i  =  %i/%i\n\n",
          a.numer,a.denom,b.numer,b.denom,
          prod.numer,prod.denom);
   make_proper(&prod,&pprod);
   printf("Proper Fraction:");
   if (abs(pprod.whole)>0)
     { printf("%i(%i/%i)",pprod.whole,pprod.numer,pprod.denom); }
   else
     { printf("%i/%i",pprod.numer,pprod.denom); }
   printf("\n");

   printf("\n %i/%i  /  %i/%i  =  %i/%i\n\n",
          a.numer,a.denom,b.numer,b.denom,
          quot.numer,quot.denom);
   make_proper(&quot,&pquot);
   printf("Proper Fraction:");
   if (abs(pquot.whole)>0)
     { printf("%i(%i/%i)",pquot.whole,pquot.numer,pquot.denom); }
   else
    { printf("%i/%i",pquot.numer,pquot.denom); }
   printf("\n");

}

void fraction_add(const fraction* a,
                  const fraction* b,
                  fraction* sum)
{
   sum->denom = a->denom * b->denom;
   sum->numer = a->numer * b->denom + b->numer * a->denom;

   void fraction_reduce(fraction* sum);
   fraction_reduce(sum);
}

void fraction_substract(const fraction* a,
                  const fraction* b,
                  fraction* diff)
{
   diff->denom = a->denom * b->denom;
   diff->numer = a->numer * b->denom - b->numer * a->denom;

   void fraction_reduce(fraction* diff);
   fraction_reduce(diff);
}

void fraction_mult(const fraction* a,
                  const fraction* b,
                  fraction* prod)
{
   prod->denom = a->denom * b->denom;
   prod->numer = a->numer * b->numer;

   void fraction_reduce(fraction* prod);
   fraction_reduce(prod);
}

void fraction_divide(const fraction* a,
                  const fraction* b,
                  fraction* quot)
{
   quot->denom = a->denom * b->numer;
   quot->numer = a->numer * b->denom;

   void fraction_reduce(fraction* quot);
   fraction_reduce(quot);
}


void make_proper(const fraction* in,
                      properfraction* out)
{
  if (in->numer*in->denom>0){
    if (in->numer>0){
      out->whole=(int)in->numer/in->denom;
      out->numer=in->numer%in->denom;
      out->denom=in->denom;
    }else{
      int Pn = - in->numer;
      int Pd = - in->denom;
      out->whole=(int)Pn/Pd;
      out->numer=Pn%Pd;
      out->denom=Pd;
    }
  } else{
    if (-(int)in->numer/in->denom>0){
      if (in->numer>0){
        int Pd = - in->denom;
        out->whole=-(int)in->numer/Pd;
        out->numer=in->numer%Pd;
        out->denom=Pd;
      }else{
        int Pn = - in->numer;
        out->whole = - (int)Pn/in->denom;
        out->numer = Pn%in->denom;
        out->denom = in->denom;
      }
    }else{
      out->whole=0;
      if (in->denom<0){
        out->numer = - in->numer;
        out->denom = - in->denom;
      }else{
        out->numer = in->numer;
        out->denom = in->denom;
      }
    }
  }

}

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
