#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Creates a deep copy of a source matrix.
 *
 * Allocates new memory and copies all values from the source matrix to the
 * destination matrix. The destination matrix must be uninitialized or properly
 * destroyed before calling this function to avoid memory leaks.
 *
 * @param src Source matrix to copy from
 * @param dst Destination matrix to copy to (will be allocated)
 *
 * @note The destination matrix's memory is allocated and must be destroyed
 *       with destroy_matrix() when no longer needed.
 * @warning If allocation fails, dst->values will be NULL and the function
 *          returns without copying any data.
 *
 * @see allocate_matrix
 * @see destroy_matrix
 */
void copy_matrix(Matrix src, Matrix *dst) {
  dst->n = src.n;
  dst->values = allocate_matrix(src.n);
  if (!dst->values)
    return;

  for (int i = 0; i < src.n; i++) {
    memcpy(dst->values[i], src.values[i], sizeof(int) * src.n);
  }
}

/**
 * @brief Allocates memory for an n×n integer matrix.
 *
 * Creates a 2D array of integers with dimensions n×n. Each row is individually
 * allocated to enable efficient cache usage in blocked algorithms.
 *
 * @param n Dimension of the square matrix (n × n)
 * @return Pointer to the allocated matrix, or NULL if allocation fails
 *
 * @note The matrix is allocated as an array of n pointers, each pointing to
 *       an array of n integers.
 * @warning If allocation fails at any row, all previously allocated rows are
 *          freed and NULL is returned.
 *
 * Example usage:
 * @code
 * int **matrix = allocate_matrix(100);
 * if (!matrix) {
 *     // Handle allocation failure
 * }
 * @endcode
 *
 * @see destroy_matrix
 */
int **allocate_matrix(int n) {
  int **matrix = (int **)malloc(n * sizeof(int *));
  if (!matrix)
    return NULL;

  for (int i = 0; i < n; i++) {
    matrix[i] = (int *)malloc(n * sizeof(int));
    if (!matrix[i]) {
      // Cleanup on failure
      for (int j = 0; j < i; j++) {
        free(matrix[j]);
      }
      free(matrix);
      return NULL;
    }
  }
  return matrix;
}

/**
 * @brief Safely deallocates a 2D integer buffer.
 *
 * Frees all memory associated with a 2D integer array allocated by
 * allocate_matrix(). This function is safe to call with NULL pointers
 * or partially allocated buffers.
 *
 * @param buf Pointer to the 2D buffer to destroy
 * @param n Number of rows in the buffer (must match allocation size)
 *
 * @note This function can handle buffers that were only partially allocated
 *       due to allocation failures.
 * @warning The parameter n must accurately reflect the number of allocated
 *          rows to avoid memory access violations.
 *
 * @see allocate_matrix
 * @see destroy_matrix
 */
void destroy_buf(int **buf, int n) {
  if (!buf)
    return;

  for (int i = 0; i < n; i++) {
    free(buf[i]);
  }
  free(buf);
}

/**
 * @brief Destroys a Matrix structure and frees its memory.
 *
 * Completely deallocates all memory associated with a Matrix structure,
 * including all row arrays. Safe to call on partially allocated matrices
 * or with NULL values.
 *
 * @param m Matrix structure to destroy
 *
 * @note This function checks for NULL pointers before attempting to free,
 *       making it safe to call multiple times or on uninitialized matrices.
 * @warning After calling this function, the Matrix structure should not be
 *          used unless reinitialized.
 *
 * Example usage:
 * @code
 * Matrix m = read_input_matrix_from_file("input.txt");
 * // Use matrix...
 * destroy_matrix(m);
 * @endcode
 *
 * @see allocate_matrix
 * @see copy_matrix
 */
void destroy_matrix(Matrix m) {
  if (m.values) {
    for (int i = 0; i < m.n; i++)
      free(m.values[i]);
    free(m.values);
  }
}

/**
 * @brief Verifies that two matrices contain identical APSP results.
 *
 * Compares two matrices element-by-element to ensure they contain the same
 * shortest path distances. Uses assertions to validate matrix dimensions
 * and all element values.
 *
 * @param a First matrix (computed result)
 * @param b Second matrix (expected result)
 * @param title Descriptive name for the verification (used in output message)
 *
 * @note If verification fails, the program will terminate with an assertion
 *       error showing the first mismatching element.
 * @note On success, prints a confirmation message with the matrix size.
 * @warning Both matrices must have the same dimensions n×n.
 *
 * Example output on success:
 * @code
 * "floyd_warshall: Assertion for matrix of size 100 was successful!"
 * @endcode
 *
 * @see Matrix
 * @assert a.n == b.n
 * @assert For all i,j: a.values[i][j] == b.values[i][j]
 */
void assert_apsp(Matrix a, Matrix b, char *title) {
  assert(a.n == b.n);
  for (int i = 0; i < a.n; i++) {
    for (int j = 0; j < a.n; j++) {
      assert(a.values[i][j] == b.values[i][j]);
    }
  }
  printf("%s: Assertion for matrix of size %d was successful!\n", title, a.n);
}
