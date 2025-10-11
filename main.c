#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int **values;
    int n;
} Matrix;

Matrix read_matrix_from_file(char* file_path);
void print_matrix(Matrix matrix);
void destroy_matrix(Matrix matrix);
int special_matrix_mul(Matrix a, Matrix b, Matrix *buf);

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Not enough arguments!\n");
        return 1;
    }
    Matrix matrix = read_matrix_from_file(argv[1]);
    print_matrix(matrix);

    Matrix matrix_mul;
    special_matrix_mul(matrix, matrix, &matrix_mul);
    print_matrix(matrix_mul);

    destroy_matrix(matrix);
    destroy_matrix(matrix_mul);
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
        }
    }

    fclose(file);

    return (Matrix) { values, n };
}

void print_matrix(Matrix matrix) {
    printf("Matrix of size (%dx%d)\n", matrix.n, matrix.n);
    for (int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.n; j++) {
            printf("%d ", matrix.values[i][j]);
        }
        printf("\n");
    }
}

void destroy_matrix(Matrix matrix) {
    if (matrix.values != NULL) free(matrix.values);
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
