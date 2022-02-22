#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fraction.h"
#include "properfraction.h"

int main()
{
   fraction a,b,sum,diff,prod,quot;

   a.numer = 1; a.denom = 10;
   b.numer = 5; b.denom = 3;


   properfraction psum,pdiff,pprod,pquot;

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
