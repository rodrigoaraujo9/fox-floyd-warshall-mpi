#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double **data;
  int size;
} Matrix;

int main() {
  // variable declarations
  int size, i, j, **matrix;

  // parse size from file
  FILE *file = fopen("matrix_examples/input1200", "r");
  if (file == NULL) {
    perror("error opening file");
    return 1;
  }
  if (fscanf(file, "%d", &size) != 1) {
    fprintf(stderr, "could not read size\n");
    fclose(file);
    return 1;
  }
  // debug
  printf("size: %d \n", size);

  // allocate mem for matrix
  matrix = malloc(size * sizeof(int *));
  for (i = 0; i < size; i++) {
    matrix[i] = malloc(size * sizeof(int));
  }

  // read matrix from file
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      if (fscanf(file, "%d", &matrix[i][j]) != 1) {
        fprintf(stderr, "could not read [%d][%d]\n", i, j);
        fclose(file);
        return 1;
      } else if (i == j && matrix[i][j] != 0) {
        fprintf(stderr,
                "distance greater than 0 between same vertice in [%d][%d]\n", i,
                j);
        return 1;
      }
    }
  }

  fclose(file);

  // debug
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      printf("%d ", matrix[i][j]);
    }
    printf("\n");
  }

  // free data from matrix
  for (int i = 0; i < size; i++) {
    free(matrix[i]);
  }
  free(matrix);

  return 0;
}
