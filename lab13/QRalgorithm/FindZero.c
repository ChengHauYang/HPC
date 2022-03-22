#include <stdlib.h>
#include <math.h>
#include "trimatrix.h"

int FindZero(trimatrix* T)
{
  const int N = T->rows;
  int k=0;
  int mfound=0;

/*
 you are having difficulty getting the QR algorithm to converge, consider changing the
     tolerance TOL by editing the FindZero.c file and changing the line
*/

//  double TOL = 1.0e-15;
  double TOL = 1.0e-10;


  while (mfound==0 && k<(N-1))
    {
      k=k+1;
      if (fabs(tget(T,k,k+1))<TOL)
        { mfound = 1; }
    }

  if (mfound==0)
    {
      k=k+1;
    }

  return k;
}
