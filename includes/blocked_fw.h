#ifndef BLOCKED_FW_H
#define BLOCKED_FW_H

#include "types.h"

int **get_block(int **matrix, int b_row, int b_col, int b, int n);

int blocked_floyd_warshall_apsp(Matrix w, Matrix *buf, int b);
int blocked_floyd_warshall_p_apsp(Matrix w, Matrix *buf, int b);

void floyd_kernel(int **C, int **A, int **B, int b);

#endif // BLOCKED_FW_H
