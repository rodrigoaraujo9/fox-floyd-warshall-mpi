#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int **allocate_matrix(int n) {
  int **matrix = (int **)malloc(n * sizeof(int *));
  if (!matrix)
    return NULL;

  for (int i = 0; i < n; i++) {
    matrix[i] = (int *)malloc(n * sizeof(int));
    if (!matrix[i]) {
      // Cleanup on failure
      for (int j = 0; j < i; j++) {
        free(matrix[j]);
      }
      free(matrix);
      return NULL;
    }
  }
  return matrix;
}

void destroy_buf(int **buf, int n) {
  if (!buf)
    return;

  for (int i = 0; i < n; i++) {
    free(buf[i]);
  }
  free(buf);
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
