#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "vector.h"

vector new_vector(const int size)
{
  vector vec;
  vec.size = size;
  assert(size>0);
  vec.val = (double*)malloc(sizeof(double)*size);

  return vec;
}

void delete_vector(vector* v)
{
  free(v->val);
}

void copy_vector(const vector* vec_in, vector* vec_out)
{
  const int N = vec_in->size;

  for (int i=0; i<N; i++)
    { vec_out->val[i] = vec_in->val[i]; }
}
