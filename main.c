#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
  int **values;
  int n;
} Matrix;

typedef struct {
  char *i;
  char *o;
} IO_Files;

IO_Files io_files[] = {
    {.i = "matrix_examples/input5", .o = "matrix_examples/output5"},
    {.i = "matrix_examples/input6", .o = "matrix_examples/output6"},
    {.i = "matrix_examples/input300", .o = "matrix_examples/output300"},
    {.i = "matrix_examples/input600", .o = "matrix_examples/output600"},
    /*
    {.i = "matrix_examples/input900", .o = "matrix_examples/output900"},
    {.i = "matrix_examples/input1200", .o = "matrix_examples/output1200"},
     */
};

Matrix read_input_matrix_from_file(char *file_path);
Matrix read_output_matrix_from_file(char *file_path, int n);

void print_matrix(Matrix matrix, char *title);
void copy_matrix(Matrix matrix_to_copy, Matrix *buf);
void destroy_matrix(Matrix matrix);

int special_matrix_mul(Matrix a, Matrix b, Matrix *buf);

int slow_apsp(Matrix w, Matrix *buf);
int repeated_squaring_apsp(Matrix w, Matrix *buf);
int floyd_warshall_apsp(Matrix w, Matrix *buf);

void assert_apsp(Matrix a, Matrix b, char *title);

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Not enough arguments!\n");
    return 1;
  }

  clock_t start, end;

  Matrix input_file_matrix = read_input_matrix_from_file(argv[1]);
  // print_matrix(matrix, "Input Matrix");

  // Matrix matrix_slow_apsp;
  // start = clock();
  // if (slow_apsp(matrix, &matrix_slow_apsp) != 0) {
  //     fprintf(stderr, "Failed to apsp!\n");
  // }
  // end = clock();
  // print_matrix(matrix_slow_apsp, "Slow APSP");
  // printf("Slow APSP: Speed = %f\n", (float)(end - start) / CLOCKS_PER_SEC);

  Matrix matrix_rs_apsp;
  start = clock();
  if (repeated_squaring_apsp(input_file_matrix, &matrix_rs_apsp) != 0) {
    fprintf(stderr, "Failed to apsp!\n");
  }
  end = clock();
  // print_matrix(matrix_rs_apsp, "Repeated Squaring APSP");
  printf("Repeated Squaring APSP: Speed = %f\n",
         (float)(end - start) / CLOCKS_PER_SEC);

  Matrix matrix_fw_apsp;
  start = clock();
  if (floyd_warshall_apsp(input_file_matrix, &matrix_fw_apsp) != 0) {
    fprintf(stderr, "Failed to apsp!\n");
  }
  end = clock();
  // print_matrix(matrix_fw_apsp, "Floyd Warshall APSP");
  printf("Floyd Warshall APSP: Speed = %f\n",
         (float)(end - start) / CLOCKS_PER_SEC);

  Matrix output_file_matrix =
      read_output_matrix_from_file(argv[2], input_file_matrix.n);

  assert_apsp(matrix_rs_apsp, output_file_matrix, "Repeated Squaring APSP");
  assert_apsp(matrix_fw_apsp, output_file_matrix, "Floyd Warshall APSP");

  destroy_matrix(input_file_matrix);
  destroy_matrix(output_file_matrix);
  // destroy_matrix(matrix_slow_apsp);
  destroy_matrix(matrix_rs_apsp);
  destroy_matrix(matrix_fw_apsp);
  return 0;
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

Matrix read_input_matrix_from_file(char *file_path) {
  FILE *file;
  int n, **values;

  // open the file
  if ((file = fopen(file_path, "r")) == NULL) {
    return (Matrix){values, n};
  }

  // scan the size of the rows and columns
  fscanf(file, "%d", &n);

  // allocate space for each column of the matrix
  values = (int **)malloc(sizeof(int *) * n);

  for (int i = 0; i < n; i++) {
    // allocate space for each line of the matrix
    values[i] = (int *)malloc(sizeof(int) * n);
    for (int j = 0; j < n; j++) {
      // read the value and assign it to the correct position in the matrix
      fscanf(file, "%d", &values[i][j]);
      if (i != j && values[i][j] == 0) {
        values[i][j] = INT_MAX;
      }
    }
  }

  fclose(file);

  return (Matrix){values, n};
}

Matrix read_output_matrix_from_file(char *file_path, int n) {
  FILE *file;
  int **values;

  // open the file
  if ((file = fopen(file_path, "r")) == NULL) {
    return (Matrix){values, n};
  }

  // allocate space for each column of the matrix
  values = (int **)malloc(sizeof(int *) * n);

  for (int i = 0; i < n; i++) {
    // allocate space for each line of the matrix
    values[i] = (int *)malloc(sizeof(int) * n);
    for (int j = 0; j < n; j++) {
      // read the value and assign it to the correct position in the matrix
      fscanf(file, "%d", &values[i][j]);
      if (i != j && values[i][j] == 0) {
        values[i][j] = INT_MAX;
      }
    }
  }

  fclose(file);

  return (Matrix){values, n};
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

void copy_matrix(Matrix matrix_to_copy, Matrix *buf) {
  (*buf).n = matrix_to_copy.n;
  (*buf).values = (int **)malloc(sizeof(int *) * matrix_to_copy.n);

  for (int i = 0; i < matrix_to_copy.n; i++) {
    (*buf).values[i] = (int *)malloc(sizeof(int) * matrix_to_copy.n);
    for (int j = 0; j < matrix_to_copy.n; j++) {
      (*buf).values[i][j] = matrix_to_copy.values[i][j];
    }
  }
}

void destroy_matrix(Matrix matrix) {
  if (matrix.values != NULL) {
    for (int i = 0; i < matrix.n; i++) {
      if (matrix.values[i] != NULL)
        free(matrix.values[i]);
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

  *buf = (Matrix){(int **)malloc(sizeof(int *) * a.n), a.n};
  int help;

  for (int i = 0; i < (*buf).n; i++) {
    (*buf).values[i] = (int *)malloc(sizeof(int) * (*buf).n);
    for (int j = 0; j < (*buf).n; j++) {
      (*buf).values[i][j] = INT_MAX;
      for (int k = 0; k < (*buf).n; k++) {
        if (a.values[i][k] != INT_MAX && b.values[k][j] != INT_MAX) {
          help = (a.values[i][k] + b.values[k][j]);
          if ((*buf).values[i][j] > help) {
            (*buf).values[i][j] = help;
          }
        }
      }
    }
  }

  return 0;
}

// 3 implementation sequential
int slow_apsp(Matrix w, Matrix *buf) {
  Matrix d_next;

  copy_matrix(w, buf);

  // from d(2) to d(n-1)
  for (int k = 2; k < w.n - 1; k++) {
    special_matrix_mul(*buf, w, &d_next);
    if (d_next.values == NULL) {
      fprintf(stderr, "special_matrix_mul failed at k=%d\n", k);
      return 1;
    }
    destroy_matrix(*buf);
    *buf = d_next;
  }
  return 0;
}

// 3.1 implementation sequential
int repeated_squaring_apsp(Matrix w, Matrix *buf) {
  int m;
  Matrix d_2m;
  copy_matrix(w, buf);

  if ((*buf).values == NULL) {
    fprintf(stderr, "Failed to copy initial matrix\n");
    return 1;
  }

  m = 1;

  while (m < w.n - 1) {
    special_matrix_mul((*buf), (*buf), &d_2m);
    if (d_2m.values == NULL) {
      destroy_matrix(*buf);
      fprintf(stderr, "special_matrix_mul failed at m=%d\n", m);
      return 1;
    }

    m *= 2;
    destroy_matrix(*buf);
    (*buf) = d_2m;
  }
  return 0;
}

int floyd_warshall_apsp(Matrix w, Matrix *buf) {
  int new_d, *ri, *rk, rik;
  copy_matrix(w, buf);
  if ((*buf).values == NULL) {
    fprintf(stderr, "Failed to copy initial matrix\n");
    return 1;
  }

  for (int k = 0; k < w.n; k++) {
    rk = (*buf).values[k];
    for (int i = 0; i < w.n; i++) {
      ri = (*buf).values[i];
      rik = ri[k];
      if (rik == INT_MAX)
        continue;
      for (int j = 0; j < w.n; j++) {
        if (rk[j] != INT_MAX) {
          new_d = rik + rk[j];
          if (ri[j] > new_d) {
            ri[j] = new_d;
          }
        }
      }
    }
  }
  return 0;
}
