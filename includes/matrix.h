#ifndef MATRIX_H
#define MATRIX_H

#include "types.h"

void copy_matrix(Matrix src, Matrix *dst);
void destroy_matrix(Matrix m);
void assert_apsp(Matrix a, Matrix b, char *title);
int allocate_matrix(int **buf, int n);

#endif // MATRIX_H
