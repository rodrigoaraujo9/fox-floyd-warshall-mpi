#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

void copy_matrix(Matrix src, Matrix *dst) {
  dst->n = src.n;
  dst->values = (int **)malloc(sizeof(int *) * src.n);
  for (int i = 0; i < src.n; i++) {
    dst->values[i] = (int *)malloc(sizeof(int) * src.n);
    for (int j = 0; j < src.n; j++) {
      dst->values[i][j] = src.values[i][j];
    }
  }
}

void destroy_matrix(Matrix m) {
  if (m.values) {
    for (int i = 0; i < m.n; i++)
      free(m.values[i]);
    free(m.values);
  }
}

void assert_apsp(Matrix a, Matrix b, char *title) {
  assert(a.n == b.n);
  for (int i = 0; i < a.n; i++) {
    for (int j = 0; j < a.n; j++) {
      assert(a.values[i][j] == b.values[i][j]);
    }
  }
  printf("%s: Assertion for matrix of size %d was successful!\n", title, a.n);
}
