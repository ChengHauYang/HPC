#ifndef __VECTOR_H__
#define __VECTOR_H__

typedef struct vector vector;
struct vector
{ int size; double* val; };

// Shortcut evaluate functions
#define vget(vec,i) vec.val[(i-1)]
#define vgetp(vec,i) vec->val[(i-1)]

// Function declarations
vector new_vector(const int size);
void delete_vector(vector* v);
void copy_vector(const vector* vec_in, vector* vec_out);

#endif
