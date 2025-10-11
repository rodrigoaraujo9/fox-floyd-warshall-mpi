#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int **values;
    int n;
} Matrix;

Matrix read_matrix_from_file(char* file_path);
void print_matrix(Matrix matrix);
void destroy_matrix(Matrix matrix);

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Not enough arguments");
        return 1;
    }
    Matrix matrix = read_matrix_from_file(argv[1]);
    print_matrix(matrix);
    destroy_matrix(matrix);
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
