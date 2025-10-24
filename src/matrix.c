#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

int allocate_matrix(int **buf, int n) {
    if (buf == NULL) free(buf);
    buf = (int**) malloc(sizeof(int*) * n);
    for (int i = 0; i < n; i++) {
        buf[i] = (int*) malloc(sizeof(int) * n);
    }
    return 0;
}

void copy_matrix(Matrix src, Matrix *dst) {
    memcpy(&src.values, &(*dst).values, sizeof(int) * src.n * src.n);
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
