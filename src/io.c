#include "../includes/io.h"
#include "../includes/types.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

Matrix read_input_matrix_from_file(char *file_path) {
  FILE *file = fopen(file_path, "r");
  if (!file)
    return (Matrix){NULL, 0};

  int n;
  if (fscanf(file, "%d", &n) != 1) {
    fclose(file);
    return (Matrix){NULL, 0};
  }

  int **values = (int **)malloc(sizeof(int *) * n);
  for (int i = 0; i < n; i++) {
    values[i] = (int *)malloc(sizeof(int) * n);
    for (int j = 0; j < n; j++) {
      if (fscanf(file, "%d", &values[i][j]) != 1) {
        fclose(file);
        return (Matrix){NULL, 0};
      }
      if (i != j && values[i][j] == 0)
        values[i][j] = INT_MAX;
    }
  }

  fclose(file);
  return (Matrix){values, n};
}

Matrix read_output_matrix_from_file(char *file_path, int n) {
  FILE *file = fopen(file_path, "r");
  if (!file)
    return (Matrix){NULL, 0};

  int **values = (int **)malloc(sizeof(int *) * n);
  for (int i = 0; i < n; i++) {
    values[i] = (int *)malloc(sizeof(int) * n);
    for (int j = 0; j < n; j++) {
      if (fscanf(file, "%d", &values[i][j]) != 1) {
        fclose(file);
        return (Matrix){NULL, 0};
      }
      if (i != j && values[i][j] == 0)
        values[i][j] = INT_MAX;
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
      } else {
        printf("%d ", matrix.values[i][j]);
      }
    }
    printf("\n");
  }
}
