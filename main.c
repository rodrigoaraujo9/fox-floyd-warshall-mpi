#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int **allocate_matrix(int size) {
  int **matrix = malloc(size * sizeof(int *));
  if (matrix == NULL) {
    return NULL;
  }
  for (int i = 0; i < size; i++) {
    matrix[i] = malloc(size * sizeof(int));
    if (matrix[i] == NULL) {
      for (int k = 0; k < i; k++) {
        free(matrix[k]);
      }
      free(matrix);
      return NULL;
    }
  }
  return matrix;
}

void free_matrix(int **matrix, int size) {
  if (matrix == NULL) {
    return;
  }
  for (int i = 0; i < size; i++) {
    free(matrix[i]);
  }
  free(matrix);
}

void print_matrix(int **matrix, int size, const char *label) {
  if (label != NULL) {
    printf("%s\n", label);
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (matrix[i][j] == INT_MAX) {
        printf("INF ");
      } else {
        printf("%d   ", matrix[i][j]);
      }
    }
    printf("\n");
  }
}

int read_size(FILE *file, int *size) {
  if (fscanf(file, "%d", size) != 1) {
    fprintf(stderr, "could not read size\n");
    return 0;
  }
  return 1;
}

int read_matrix(FILE *file, int **matrix, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (fscanf(file, "%d", &matrix[i][j]) != 1) {
        fprintf(stderr, "could not read [%d][%d]\n", i, j);
        return 0;
      }
      // Validate diagonal elements
      if (i == j && matrix[i][j] != 0) {
        fprintf(stderr,
                "distance greater than 0 between same vertice in [%d][%d]\n", i,
                j);
        return 0;
      }
      if (i != j && matrix[i][j] == 0) {
        matrix[i][j] = INT_MAX;
      }
    }
  }
  return 1;
}

int **parse_matrix_file(const char *filename, int *out_size) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    perror("error opening file");
    return NULL;
  }

  int size;
  if (!read_size(file, &size)) {
    fclose(file);
    return NULL;
  }

  int **matrix = allocate_matrix(size);
  if (matrix == NULL) {
    fprintf(stderr, "memory allocation failed\n");
    fclose(file);
    return NULL;
  }

  if (!read_matrix(file, matrix, size)) {
    free_matrix(matrix, size);
    fclose(file);
    return NULL;
  }

  fclose(file);
  *out_size = size;
  return matrix;
}

// given algorythm for matrix mul using dynamic programming
int **special_matrix_mul(int **a, int a_size, int **b, int b_size) {
  if (a_size != b_size) {
    fprintf(stderr, "matrix sizes must be same: %d != %d\n", a_size, b_size);
    return NULL;
  }

  int **c, i, j, k, help, n = a_size;

  c = allocate_matrix(a_size);
  if (c == NULL) {
    fprintf(stderr, "mem alloc failed in special_matrix_mul\n");
    return NULL;
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      c[i][j] = INT_MAX;
      for (k = 0; k < n; k++) {
        if (a[i][k] != INT_MAX && b[k][j] != INT_MAX) {
          help = (a[i][k] + b[k][j]);
          if (c[i][j] > help) {
            c[i][j] = help;
          }
        }
      }
    }
  }
  return c;
}

int **copy_matrix(int **m, int size) {
  int **copy = allocate_matrix(size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      copy[i][j] = m[i][j];
    }
  }
  return copy;
}

int **slow_apsp(int **w, int size) {
  int k, **d_next;

  int **d = copy_matrix(w, size);

  // from d(2) to d(n-1)
  for (k = 2; k < size; k++) {
    d_next = special_matrix_mul(d, size, w, size);
    if (d_next == NULL) {
      free_matrix(d, size);
      fprintf(stderr, "special_matrix_mul failed at k=%d\n", k);
      return NULL;
    }
    free_matrix(d, size);
    d = d_next;
  }
  return d;
}

int main() {
  int size;
  int **w = parse_matrix_file("matrix_examples/input5", &size);
  if (w == NULL) {
    return 1;
  }

  printf("size of matrix: %d\n", size);

  print_matrix(w, size, "w");
  int **d = slow_apsp(w, size);
  if (d == NULL) {
    fprintf(stderr, "apsp failed\n");
    free_matrix(w, size);
    return 1;
  }

  print_matrix(d, size, "final d");

  free_matrix(w, size);
  free_matrix(d, size);
  return 0;
}
