#include "../includes/io.h"
#include "../includes/types.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Reads a weighted graph adjacency matrix from a file
 *
 * Reads a matrix from a text file where the first line contains the
 * dimension n, followed by n lines of n space-separated integers
 * representing edge weights. The function automatically converts
 * 0 values to INT_MAX (infinity) to represent missing edges, except
 * on the diagonal where 0 represents self-loops.
 *
 * File format example for a 3×3 matrix:
 * @code
 * 3
 * 0 5 0
 * 0 0 3
 * 2 0 0
 * @endcode
 *
 * This represents:
 * - Edge 0→1 with weight 5
 * - Edge 1→2 with weight 3
 * - Edge 2→0 with weight 2
 * - All other edges are missing (infinity)
 * - Self-loops have weight 0
 *
 * @param file_path Path to the input file
 *
 * @return Matrix structure containing the loaded matrix, or {NULL, 0} on failure
 *
 * @note Returns {NULL, 0} if file cannot be opened or format is invalid
 * @note Off-diagonal zeros are converted to INT_MAX (no edge)
 * @note Diagonal entries remain 0 (self-loops)
 * @note Caller must free the returned matrix using destroy_matrix()
 *
 * @warning Memory is allocated and must be freed by caller
 * @warning Check matrix.values != NULL before using returned matrix
 *
 * Error handling:
 * - File not found: returns {NULL, 0}
 * - Invalid format: returns {NULL, 0} and closes file
 * - Insufficient data: returns {NULL, 0} and closes file
 *
 * @see read_output_matrix_from_file for reading pre-computed results
 * @see destroy_matrix for cleanup
 */
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

/**
 * @brief Reads a pre-computed output matrix from a file
 *
 * Reads a matrix of shortest path distances from a text file. Similar to
 * read_input_matrix_from_file but requires the dimension n as a parameter
 * instead of reading it from the file. This is used for loading expected
 * results for verification.
 *
 * The file contains n lines of n space-separated integers. Like the input
 * format, 0 values are converted to INT_MAX except on the diagonal.
 *
 * @param file_path Path to the output file
 * @param n Matrix dimension (must be known in advance)
 *
 * @return Matrix structure containing the loaded matrix, or {NULL, 0} on failure
 *
 * @note Unlike read_input_matrix_from_file, dimension n is a parameter
 * @note Off-diagonal zeros are converted to INT_MAX (no path)
 * @note Diagonal entries remain 0 (distance to self)
 * @note Caller must free the returned matrix using destroy_matrix()
 *
 * @warning Memory is allocated and must be freed by caller
 * @warning Check matrix.values != NULL before using returned matrix
 *
 * Use case:
 * This function is typically used for loading golden reference outputs
 * to verify the correctness of APSP algorithm implementations.
 *
 * @see read_input_matrix_from_file for reading input graphs
 * @see destroy_matrix for cleanup
 */
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

/**
 * @brief Prints a matrix to standard output
 *
 * Displays a matrix with a title header showing the dimensions. For display
 * purposes, INT_MAX values (infinity) are printed as 0 to match the input
 * file format convention.
 *
 * Output format:
 * @code
 * Title (3x3)
 * 0 5 0
 * 0 0 3
 * 2 0 0
 * @endcode
 *
 * @param matrix The matrix to print
 * @param title Descriptive title to display above the matrix
 *
 * @note INT_MAX values are displayed as 0 (no path)
 * @note Format matches the input file format for consistency
 * @note Useful for debugging and verifying results
 *
 * Usage example:
 * @code
 * Matrix result;
 * floyd_warshall_apsp(input, &result);
 * print_matrix(result, "Shortest Paths");
 * @endcode
 *
 * @see read_input_matrix_from_file for the inverse operation
 */
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

/**
 * @brief Flushes CPU cache for accurate performance measurements
 *
 * Allocates and accesses a large block of memory (20MB) to evict cached
 * data from the CPU caches. This ensures that performance measurements
 * start with a cold cache, providing more consistent and reproducible
 * timing results across multiple benchmark runs.
 *
 * The function:
 * 1. Allocates 20MB of memory
 * 2. Writes to every byte to ensure it's loaded into cache
 * 3. Reads the last byte to prevent compiler optimization
 * 4. Frees the memory
 *
 * @note 20MB is typically larger than L3 cache on most modern CPUs
 * @note Uses volatile to prevent compiler from optimizing away the access
 * @note Silently fails if allocation fails (benchmarking continues)
 *
 * @warning Should only be used in benchmarking code, not production
 * @warning Adds overhead to measurements (call before timing starts)
 *
 * Usage in benchmarks:
 * @code
 * for (int trial = 0; trial < num_trials; trial++) {
 *     flush_cache();  // Start with cold cache
 *     start_timer();
 *     algorithm();
 *     stop_timer();
 * }
 * @endcode
 *
 * Cache sizes (typical):
 * - L1: 32-64 KB per core
 * - L2: 256 KB - 1 MB per core
 * - L3: 8-32 MB shared
 *
 * Purpose:
 * Without cache flushing, the first run may be slower than subsequent
 * runs due to cache warming effects, leading to inconsistent measurements.
 *
 * @see Performance benchmarking best practices
 */
void flush_cache() {
  const size_t size = 20 * 1024 * 1024; // 20 MB
  char *c = malloc(size);
  if (c) {
    for (size_t i = 0; i < size; i++) {
      c[i] = i;
    }
    volatile char dummy = c[size - 1];
    free(c);
    (void)dummy;
  }
}
