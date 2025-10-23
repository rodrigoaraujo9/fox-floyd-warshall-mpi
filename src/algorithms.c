#include "../includes/algorithms.h"
#include "../includes/matrix.h"
#include "../includes/types.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

int special_matrix_mul(Matrix a, Matrix b, Matrix *buf) {
  if (a.n != b.n) {
    fprintf(stderr, "matrix sizes must be same: %d != %d\n", a.n, b.n);
    return 1;
  }

  *buf = (Matrix){(int **)malloc(sizeof(int *) * a.n), a.n};
  for (int i = 0; i < buf->n; i++) {
    buf->values[i] = (int *)malloc(sizeof(int) * buf->n);
    for (int j = 0; j < buf->n; j++) {
      buf->values[i][j] = INT_MAX;
      for (int k = 0; k < buf->n; k++) {
        if (a.values[i][k] != INT_MAX && b.values[k][j] != INT_MAX) {
          int sum = a.values[i][k] + b.values[k][j];
          if (buf->values[i][j] > sum)
            buf->values[i][j] = sum;
        }
      }
    }
  }
  return 0;
}

int slow_apsp(Matrix w, Matrix *buf) {
  Matrix d_next;
  copy_matrix(w, buf);
  for (int k = 2; k < w.n - 1; k++) {
    if (special_matrix_mul(*buf, w, &d_next) != 0)
      return 1;
    destroy_matrix(*buf);
    *buf = d_next;
  }
  return 0;
}

int repeated_squaring_apsp(Matrix w, Matrix *buf) {
  copy_matrix(w, buf);
  if (!buf->values)
    return 1;

  int m = 1;
  while (m < w.n - 1) {
    Matrix d_2m;
    if (special_matrix_mul(*buf, *buf, &d_2m) != 0) {
      destroy_matrix(*buf);
      return 1;
    }
    m *= 2;
    destroy_matrix(*buf);
    *buf = d_2m;
  }
  return 0;
}

int floyd_warshall_apsp(Matrix w, Matrix *buf) {
  copy_matrix(w, buf);
  if (!buf->values)
    return 1;

  for (int k = 0; k < w.n; k++) {
    int *rk = buf->values[k];
    for (int i = 0; i < w.n; i++) {
      int *ri = buf->values[i];
      int rik = ri[k];
      if (rik == INT_MAX)
        continue;
      for (int j = 0; j < w.n; j++) {
        if (rk[j] != INT_MAX) {
          int nd = rik + rk[j];
          if (ri[j] > nd)
            ri[j] = nd;
        }
      }
    }
  }
  return 0;
}
