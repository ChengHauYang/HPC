#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fraction.h"

void fraction_mult(const fraction* a,
                  const fraction* b,
                  fraction* prod)
{
   prod->denom = a->denom * b->denom;
   prod->numer = a->numer * b->numer;

   void fraction_reduce(fraction* prod);
   fraction_reduce(prod);
}
