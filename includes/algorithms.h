#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "types.h"

int special_matrix_mul(Matrix a, Matrix b, Matrix *buf);
int slow_apsp(Matrix w, Matrix *buf);
int repeated_squaring_apsp(Matrix w, Matrix *buf);
int floyd_warshall_apsp(Matrix w, Matrix *buf);

#endif // ALGORITHMS_H
