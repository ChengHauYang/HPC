#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fraction.h"

void fraction_divide(const fraction* a,
                  const fraction* b,
                  fraction* quot)
{
   quot->denom = a->denom * b->numer;
   quot->numer = a->numer * b->denom;

   void fraction_reduce(fraction* quot);
   fraction_reduce(quot);
}
