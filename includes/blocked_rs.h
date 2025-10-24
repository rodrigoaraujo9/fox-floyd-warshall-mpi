#ifndef BLOCKED_RS_H
#define BLOCKED_RS_H

#include "types.h"

int blocked_repeated_squaring_apsp(Matrix w, Matrix *buf, int block_size);
int blocked_matrix_mul_kernel(int **a, int **b, int **buf, int size);
int blocked_matrix_mul(Matrix a, Matrix b, Matrix *buf, int block_size);

#endif // BLOCKED_RS_H
