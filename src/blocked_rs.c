#include "../includes/blocked_fw.h"
#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int blocked_matrix_mul_kernel(int **a, int **b, int **buf, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                if (a[i][k] != INT_MAX && b[k][j] != INT_MAX) {
                    int sum = a[i][k] + b[k][j];
                    if (buf[i][j] > sum)
                        buf[i][j] = sum;
                }
            }
        }
    }
    return 0;
}

int blocked_matrix_mul(Matrix a, Matrix b, Matrix *buf, int block_size) {
    int n = a.n;
    if (n % block_size != 0) {
      fprintf(stderr, "Block size must divide matrix size %d\n", n);
      return 1;
    }

    buf->n = n;
    buf->values = (int **)malloc(sizeof(int *) * n);
    for (int i = 0; i < n; i++) {
      buf->values[i] = (int *)malloc(sizeof(int) * n);
      for (int j = 0; j < n; j++) {
        buf->values[i][j] = INT_MAX;
      }
    }

    int B = n / block_size;

    for (int k = 0; k < B; k++) {
      for (int i = 0; i < B; i++) {
        for (int j = 0; j < B; j++) {
          int **a_block = get_block(a.values, i, k, block_size, n);
          int **b_block = get_block(b.values, k, j, block_size, n);
          int **result_block = get_block(buf->values, i, j, block_size, n);

          blocked_matrix_mul_kernel(a_block, b_block, result_block, block_size);

          free(a_block);
          free(b_block);
          free(result_block);
        }
      }
    }

    return 0;
}

int blocked_repeated_squaring_apsp(Matrix w, Matrix *buf, int block_size) {
  copy_matrix(w, buf);
  if (!buf->values)
    return 1;

  int m = 1;
  while (m < w.n - 1) {
    Matrix d_2m;
    if (blocked_matrix_mul(*buf, *buf, &d_2m, block_size) != 0) {
      destroy_matrix(*buf);
      return 1;
    }
    m *= 2;
    destroy_matrix(*buf);
    *buf = d_2m;
  }
  return 0;
}
