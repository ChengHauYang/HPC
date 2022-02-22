#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fraction.h"

void fraction_substract(const fraction* a,
                  const fraction* b,
                  fraction* diff)
{
   diff->denom = a->denom * b->denom;
   diff->numer = a->numer * b->denom - b->numer * a->denom;

   void fraction_reduce(fraction* diff);
   fraction_reduce(diff);
}
