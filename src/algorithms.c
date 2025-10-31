#include "../includes/algorithms.h"
#include "../includes/matrix.h"
#include "../includes/types.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Performs special matrix multiplication for shortest path computation
 *
 * This function implements a modified matrix multiplication where:
 * - Standard multiplication is replaced by addition
 * - Standard addition is replaced by minimum operation
 * This operation is used in dynamic programming shortest path algorithms.
 *
 * @param a First input matrix (must be square)
 * @param b Second input matrix (must be square)
 * @param buf Pointer to output matrix where result will be stored
 *
 * @return 0 on success, 1 on failure
 *
 * @note Matrix dimensions must match (a.n == b.n)
 * @note INT_MAX represents infinity (no path exists)
 * @note Memory is allocated for buf and must be freed by caller
 *
 * @warning Caller is responsible for freeing allocated memory in buf
 *
 * Time complexity: O(n³) where n is the matrix dimension
 * Space complexity: O(n²) for the output matrix
 */
int special_matrix_mul(Matrix a, Matrix b, Matrix *buf) {
  if (a.n != b.n) {
    fprintf(stderr, "matrix sizes must be same: %d != %d\n", a.n, b.n);
    return 1;
  }

  *buf = (Matrix){(int **)malloc(sizeof(int *) * a.n), a.n};
  for (int i = 0; i < buf->n; i++) {
    buf->values[i] = (int *)malloc(sizeof(int) * buf->n);
    for (int j = 0; j < buf->n; j++) {
      buf->values[i][j] = INT_MAX;
      for (int k = 0; k < buf->n; k++) {
        if (a.values[i][k] != INT_MAX && b.values[k][j] != INT_MAX) {
          int sum = a.values[i][k] + b.values[k][j];
          if (buf->values[i][j] > sum)
            buf->values[i][j] = sum;
        }
      }
    }
  }
  return 0;
}

/**
 * @brief Solves APSP using slow dynamic programming approach
 *
 * Computes all-pairs shortest paths by iteratively extending paths.
 * This algorithm computes shortest paths with at most m edges for
 * m = 2, 3, ..., n-1, where n is the number of vertices.
 *
 * @param w Input weighted adjacency matrix representing the graph
 * @param buf Pointer to output matrix where shortest paths will be stored
 *
 * @return 0 on success, 1 on failure
 *
 * @note w.values[i][j] should contain the edge weight from vertex i to j
 * @note INT_MAX represents no direct edge between vertices
 * @note Result stored in buf must be freed by caller
 *
 * Time complexity: O(n⁴) where n is the number of vertices
 * Space complexity: O(n²)
 *
 * @see special_matrix_mul
 * @see repeated_squaring_apsp for a faster alternative
 */
int slow_apsp(Matrix w, Matrix *buf) {
  Matrix d_next;
  copy_matrix(w, buf);
  for (int k = 2; k < w.n - 1; k++) {
    if (special_matrix_mul(*buf, w, &d_next) != 0)
      return 1;
    destroy_matrix(*buf);
    *buf = d_next;
  }
  return 0;
}

/**
 * @brief Solves APSP using repeated squaring optimization
 *
 * Computes all-pairs shortest paths using matrix exponentiation by squaring.
 * This is an optimized version of the dynamic programming approach that
 * reduces the number of matrix multiplications from O(n) to O(log n).
 *
 * @param w Input weighted adjacency matrix representing the graph
 * @param buf Pointer to output matrix where shortest paths will be stored
 *
 * @return 0 on success, 1 on failure
 *
 * @note This algorithm squares the distance matrix repeatedly
 * @note Computes shortest paths with at most 2^m edges in iteration m
 * @note Result stored in buf must be freed by caller
 *
 * Time complexity: O(n³ log n) where n is the number of vertices
 * Space complexity: O(n²)
 *
 * @see special_matrix_mul
 * @see slow_apsp for the basic version
 * @see floyd_warshall_apsp for an even faster algorithm
 */
int repeated_squaring_apsp(Matrix w, Matrix *buf) {
  copy_matrix(w, buf);
  if (!buf->values)
    return 1;

  int m = 1;
  while (m < w.n - 1) {
    Matrix d_2m;
    if (special_matrix_mul(*buf, *buf, &d_2m) != 0) {
      destroy_matrix(*buf);
      return 1;
    }
    m *= 2;
    destroy_matrix(*buf);
    *buf = d_2m;
  }
  return 0;
}

/**
 * @brief Solves APSP using Floyd-Warshall algorithm
 *
 * Computes all-pairs shortest paths using dynamic programming with
 * intermediate vertices. This is typically the most efficient algorithm
 * for dense graphs and is simpler than the matrix multiplication approaches.
 *
 * The algorithm considers all possible intermediate vertices k and updates
 * the shortest path from i to j if the path through k is shorter.
 *
 * @param w Input weighted adjacency matrix representing the graph
 * @param buf Pointer to output matrix where shortest paths will be stored
 *
 * @return 0 on success, 1 on failure
 *
 * @note This algorithm modifies the matrix in-place for efficiency
 * @note Can handle negative edge weights (but not negative cycles)
 * @note Result stored in buf must be freed by caller
 *
 * Time complexity: O(n³) where n is the number of vertices
 * Space complexity: O(n²)
 *
 * @warning Does not detect negative weight cycles
 *
 * Algorithm overview:
 * For each intermediate vertex k:
 *   For each pair of vertices (i, j):
 *     If path i -> k -> j is shorter than current i -> j:
 *       Update shortest path i -> j
 *
 * @see repeated_squaring_apsp for an alternative approach
 */
int floyd_warshall_apsp(Matrix w, Matrix *buf) {
  copy_matrix(w, buf);
  if (!buf->values)
    return 1;

  for (int k = 0; k < w.n; k++) {
    int *rk = buf->values[k];
    for (int i = 0; i < w.n; i++) {
      int *ri = buf->values[i];
      int rik = ri[k];
      if (rik == INT_MAX)
        continue;
      for (int j = 0; j < w.n; j++) {
        if (rk[j] != INT_MAX) {
          int nd = rik + rk[j];
          if (ri[j] > nd)
            ri[j] = nd;
        }
      }
    }
  }
  return 0;
}
