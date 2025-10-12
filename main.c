#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int **values;
    int n;
} Matrix;

Matrix read_matrix_from_file(char* file_path);
void print_matrix(Matrix matrix, char *title);
void destroy_matrix(Matrix matrix);
int special_matrix_mul(Matrix a, Matrix b, Matrix *buf);
void copy_matrix(Matrix matrix_to_copy, Matrix* buf);
int slow_apsp(Matrix w, Matrix *buf);
int repeated_squaring_apsp(Matrix w, Matrix *buf);

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Not enough arguments!\n");
        return 1;
    }
    Matrix matrix = read_matrix_from_file(argv[1]);
    print_matrix(matrix, "Input Matrix");

    Matrix matrix_slow_apsp;
    if (slow_apsp(matrix, &matrix_slow_apsp) != 0) {
        fprintf(stderr, "Failed to apsp!\n");
    }
    print_matrix(matrix_slow_apsp, "Slow APSP");

    Matrix matrix_rs_apsp;
    if (slow_apsp(matrix, &matrix_rs_apsp) != 0) {
        fprintf(stderr, "Failed to apsp!\n");
    }
    print_matrix(matrix_rs_apsp, "Repeated Squaring APSP");

    destroy_matrix(matrix);
    destroy_matrix(matrix_slow_apsp);
    destroy_matrix(matrix_rs_apsp);
    return 0;
}

Matrix read_matrix_from_file(char* file_path) {
    FILE *file;
    int n, **values;

    // open the file
    if ((file = fopen(file_path, "r")) == NULL) {
        return (Matrix) {values, n};
    }

    // scan the size of the rows and columns
    fscanf(file, "%d", &n);

    // allocate space for each column of the matrix
    values = (int **) malloc(sizeof(int*) * n);

    for (int i = 0; i < n; i++) {
        // allocate space for each line of the matrix
        values[i] = (int*) malloc(sizeof(int) * n);
        for (int j = 0; j < n; j++) {
            // read the value and assign it to the correct position in the matrix
            fscanf(file, "%d", &values[i][j]);
            if (i != j && values[i][j] == 0) {
                 values[i][j] = INT_MAX;
            }
        }
    }

    fclose(file);

    return (Matrix) { values, n };
}

void print_matrix(Matrix matrix, char *title) {
    printf("%s (%dx%d)\n", title, matrix.n, matrix.n);
    for (int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.n; j++) {
            if (matrix.values[i][j] == INT_MAX) {
                printf("0 ");
                continue;
            }
            printf("%d ", matrix.values[i][j]);
        }
        printf("\n");
    }
}

void copy_matrix(Matrix matrix_to_copy, Matrix* buf) {
    (*buf).n = matrix_to_copy.n;
    (*buf).values = (int**) malloc(sizeof(int*) * matrix_to_copy.n);

    for (int i = 0; i < matrix_to_copy.n; i++) {
        (*buf).values[i] = (int*) malloc(sizeof(int) * matrix_to_copy.n);
        for (int j = 0; j < matrix_to_copy.n; j++) {
            (*buf).values[i][j] = matrix_to_copy.values[i][j];
        }
    }
}

void destroy_matrix(Matrix matrix) {
    if (matrix.values != NULL) {
        for (int i = 0; i < matrix.n; i++) {
            if (matrix.values[i] != NULL) free(matrix.values[i]);
        }
        free(matrix.values);
    }
}

// given algorythm for matrix mul using dynamic programming
int special_matrix_mul(Matrix a, Matrix b, Matrix *buf) {
  if (a.n != b.n) {
    fprintf(stderr, "matrix sizes must be same: %d != %d\n", a.n, b.n);
    return 1;
  }

  Matrix c = (Matrix) {(int**) malloc(sizeof(int*) * a.n), a.n};
  int help;

  for (int i = 0; i < c.n; i++) {
      c.values[i] = (int*) malloc(sizeof(int) * c.n);
      for (int j = 0; j < c.n; j++) {
          c.values[i][j] = INT_MAX;
          for (int k = 0; k < c.n; k++) {
              if (a.values[i][k] != INT_MAX && b.values[k][j] != INT_MAX) {
                  help = (a.values[i][k] + b.values[k][j]);
                  if (c.values[i][j] > help) {
                      c.values[i][j] = help;
                  }
              }
          }
      }
  }

  *buf = c;

  return 0;
}

// 3 implementation sequential
int slow_apsp(Matrix w, Matrix *buf) {
  Matrix d, d_next;

  copy_matrix(w, &d);

  // from d(2) to d(n-1)
  for (int k = 2; k < w.n - 1; k++) {
    special_matrix_mul(d, w, &d_next);
    if (d_next.values == NULL) {
      fprintf(stderr, "special_matrix_mul failed at k=%d\n", k);
      return 1;
    }
    destroy_matrix(d);
    d = d_next;
  }
  *buf = d;
  return 0;
}


// 3.1 implementation sequential
int repeated_squaring_apsp(Matrix w, Matrix *buf) {
  int m;
  Matrix d_2m, d_m;
  copy_matrix(w, &d_m);

  if (d_m.values == NULL) {
    fprintf(stderr, "Failed to copy initial matrix\n");
    return 1;
  }

  m = 1;

  while (m < w.n - 1) {
    special_matrix_mul(d_m, d_m, &d_2m);
    if (d_2m.values == NULL) {
      destroy_matrix(d_m);
      fprintf(stderr, "special_matrix_mul failed at m=%d\n", m);
      return 1;
    }

    m *= 2;
    destroy_matrix(d_m);
    d_m = d_2m;
  }
  *buf = d_m;
  return 0;
}
