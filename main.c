#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  char *i;
  char *o;
} i_o;

i_o match_i_o[] = {
    {.i = "matrix_examples/input5", .o = "matrix_examples/output5"},
    {.i = "matrix_examples/input6", .o = "matrix_examples/output6"},
    {.i = "matrix_examples/input300", .o = "matrix_examples/output300"},
    {.i = "matrix_examples/input600", .o = "matrix_examples/output600"},
    /*
    {.i = "matrix_examples/input900", .o = "matrix_examples/output900"},
    {.i = "matrix_examples/input1200", .o = "matrix_examples/output1200"},
     */
};

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

int **parse_output_file(const char *filename, int size) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    perror("error opening file");
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
  return matrix;
}

int **copy_matrix(int **m, int size) {
  int **copy = allocate_matrix(size);
  if (copy == NULL) {
    return NULL;
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      copy[i][j] = m[i][j];
    }
  }
  return copy;
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

// 3 implementation sequential
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

// 3.1 implementation sequential
int **repeated_squaring_apsp(int **w, int size) {
  int m, **d_2m;
  int **d_m = copy_matrix(w, size);

  if (d_m == NULL) {
    fprintf(stderr, "Failed to copy initial matrix\n");
    return NULL;
  }

  m = 1;

  while (m < size - 1) {
    d_2m = special_matrix_mul(d_m, size, d_m, size);
    if (d_2m == NULL) {
      free_matrix(d_m, size);
      fprintf(stderr, "special_matrix_mul failed at m=%d\n", m);
      return NULL;
    }

    m *= 2;
    free_matrix(d_m, size);
    d_m = d_2m;
  }
  return d_m;
}

int main() {
  int size, i, j, p, match;

  int n_pairs = sizeof(match_i_o) / sizeof(match_i_o[0]);

  for (p = 0; p < n_pairs; p++) {
    printf("making calculations for %s", match_i_o[p].i);

    int **w = parse_matrix_file(match_i_o[p].i, &size);
    int **out = parse_output_file(match_i_o[p].o, size);
    if (w == NULL || out == NULL) {
      fprintf(stderr, "failed to load pair %d\n", p);
      if (w != NULL)
        free_matrix(w, size);
      if (out != NULL)
        free_matrix(out, size);
      continue;
    }

    // printf("size of matrix: %d\n", size);
    // print_matrix(w, size, "w");

    int **d = repeated_squaring_apsp(w, size);
    if (d == NULL) {
      fprintf(stderr, "apsp failed\n");
      free_matrix(w, size);
      continue;
    }
    // print_matrix(d, size, "output computed");
    // print_matrix(out, size, "output given");

    match = 1;
    for (i = 0; i < size; i++) {
      for (j = 0; j < size; j++) {
        if (d[i][j] != out[i][j]) {
          match = 0;
          break;
        }
      }
      if (!match)
        break;
    }

    if (match) {
      printf("%s: output matches!\n", match_i_o[p].i);
    } else {
      printf("%s: output didn't match!\n", match_i_o[p].i);
    }

    free_matrix(w, size);
    free_matrix(d, size);
    free_matrix(out, size);
  }

  return 0;
}
