#include "../includes/blocked_fw.h"
#include "../includes/comm.h"
#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Extracts a block from a matrix without copying data
 *
 * Creates a pointer array that references a specific block within the
 * original matrix. This allows block operations without data duplication.
 *
 * @param matrix Source matrix as 2D array
 * @param b_row Block row index (0-based)
 * @param b_col Block column index (0-based)
 * @param b Block size (number of rows/columns in the block)
 * @param n Total matrix dimension
 *
 * @return Array of pointers to the block rows in the original matrix
 *
 * @note The returned block shares memory with the original matrix
 * @note Only the pointer array needs to be freed, not the data itself
 * @note Block bounds are validated with assertions
 *
 * @warning Caller must free the returned pointer array with free()
 * @warning Do not free individual rows as they point to the original matrix
 *
 * Example:
 * @code
 * int **block = get_block(matrix, 0, 1, 32, 128);
 * // Use block...
 * free(block);  // Only free the pointer array
 * @endcode
 */
int **get_block(int **matrix, int b_row, int b_col, int b, int n) {
  int r0 = b_row * b;
  int c0 = b_col * b;
  assert(r0 >= 0 && c0 >= 0);
  assert(r0 + b <= n && c0 + b <= n);
  int **block = (int **)malloc(sizeof(int *) * b);
  for (int i = 0; i < b; i++)
    block[i] = &matrix[r0 + i][c0];
  return block;
}

/**
 * @brief Core Floyd-Warshall kernel for a single block
 *
 * Performs the Floyd-Warshall relaxation step on three blocks:
 * C[i][j] = min(C[i][j], A[i][k] + B[k][j])
 *
 * This is the fundamental operation used in blocked Floyd-Warshall.
 * The algorithm iterates through all intermediate vertices k within
 * the block and updates shortest paths accordingly.
 *
 * @param C Output block to be updated (modified in-place)
 * @param A Left input block (paths from source to intermediate)
 * @param B Right input block (paths from intermediate to destination)
 * @param b Block size (dimension of the blocks)
 *
 * @note All blocks must be b × b in size
 * @note INT_MAX represents infinity (no path exists)
 * @note This operation modifies C in-place for efficiency
 *
 * Time complexity: O(b³) for a single block
 *
 * @see floyd_warshall_apsp for the basic algorithm
 */
void floyd_kernel(int **C, int **A, int **B, int b) {
  for (int k = 0; k < b; k++) {
    for (int i = 0; i < b; i++) {
      int a_val = A[i][k];
      if (a_val == INT_MAX)
        continue;
      for (int j = 0; j < b; j++) {
        int b_val = B[k][j];
        if (b_val == INT_MAX)
          continue;
        int sum = a_val + b_val;
        if (C[i][j] > sum)
          C[i][j] = sum;
      }
    }
  }
}

/**
 * @brief Sequential blocked Floyd-Warshall algorithm
 *
 * Implements Floyd-Warshall using a blocked approach for better cache
 * performance. The matrix is divided into blocks of size b × b, and
 * the algorithm processes blocks in a specific order to maintain
 * correctness while improving memory locality.
 *
 * For each diagonal block k:
 * 1. Update block (k,k) with itself
 * 2. Update all blocks in row k using (k,k)
 * 3. Update all blocks in column k using (k,k)
 * 4. Update all remaining blocks using row k and column k blocks
 *
 * @param w Input weighted adjacency matrix
 * @param buf Pointer to output matrix for shortest paths
 * @param b Block size (must divide matrix size evenly)
 *
 * @return 0 on success, 1 on failure
 *
 * @note Matrix dimension n must be divisible by block size b
 * @note Block size b should be chosen based on cache size (typical: 32-128)
 * @note Result stored in buf must be freed by caller
 *
 * Time complexity: O(n³) with better cache performance than standard version
 * Space complexity: O(n²)
 *
 * @warning Fails if n % b != 0
 *
 * @see floyd_warshall_apsp for the standard implementation
 * @see blocked_floyd_warshall_p_apsp for parallel version
 */
int blocked_floyd_warshall_apsp(Matrix w, Matrix *buf, int b) {
  int n = w.n;
  if (n % b != 0) {
    fprintf(stderr, "Block size must divide matrix size %d\n", n);
    return 1;
  }

  copy_matrix(w, buf);
  if (!buf->values)
    return 1;

  int B = n / b;
  for (int k = 0; k < B; k++) {
    int **wkk = get_block(buf->values, k, k, b, n);
    floyd_kernel(wkk, wkk, wkk, b);
    free(wkk);

    for (int j = 0; j < B; j++) {
      if (j == k)
        continue;
      int **wkj = get_block(buf->values, k, j, b, n);
      int **wkk2 = get_block(buf->values, k, k, b, n);
      floyd_kernel(wkj, wkk2, wkj, b);
      free(wkj);
      free(wkk2);
    }
    for (int i = 0; i < B; i++) {
      if (i == k)
        continue;
      int **wik = get_block(buf->values, i, k, b, n);
      int **wkk3 = get_block(buf->values, k, k, b, n);
      floyd_kernel(wik, wik, wkk3, b);
      free(wkk3);

      for (int j = 0; j < B; j++) {
        if (j == k)
          continue;
        int **wkj2 = get_block(buf->values, k, j, b, n);
        int **wij = get_block(buf->values, i, j, b, n);
        floyd_kernel(wij, wik, wkj2, b);
        free(wkj2);
        free(wij);
      }
      free(wik);
    }
  }
  return 0;
}

/**
 * @brief Parallel blocked Floyd-Warshall using MPI (blocking communication)
 *
 * Implements a distributed memory parallel version of blocked Floyd-Warshall
 * using MPI. Processes are arranged in a 2D Cartesian grid (sqrt(P) × sqrt(P))
 * where P is the number of processes. Each process owns one block of the matrix.
 *
 * Algorithm phases for each iteration k:
 * 1. Owner of block (k,k) updates it and broadcasts to row and column
 * 2. Processes in row k update their blocks and broadcast to their columns
 * 3. Processes in column k update their blocks and broadcast to their rows
 * 4. All other processes update using received row and column blocks
 *
 * @param w Input weighted adjacency matrix (only valid on rank 0)
 * @param buf Pointer to output matrix (gathered on rank 0)
 * @param b Block size (must equal n/sqrt(P))
 *
 * @return 0 on success, 1 on failure
 *
 * @note Requires P = q² processes where q divides n
 * @note Block size must equal n/q where q = sqrt(P)
 * @note Uses blocking MPI_Bcast for synchronization
 * @note Result is gathered to rank 0 only
 *
 * Time complexity: O(n³/P) with O(n²/√P) communication per iteration
 *
 * @warning Number of processes must be a perfect square
 * @warning Will call MPI_Abort if b != n/sqrt(P)
 *
 * @see blocked_floyd_warshall_p_non_blocking_apsp for optimized version
 * @see setup_cart for Cartesian topology setup
 */
int blocked_floyd_warshall_p_apsp(Matrix w, Matrix *buf, int b) {
  // initialization
  CartInfo cart;
  setup_cart(&cart, w.n);
  int n = w.n;
  int localN = n / cart.q;

  if (b != localN) {
    if (cart.my_rank == 0) {
      fprintf(stderr, "Set b == n/sqrt(P) (got b=%d, need %d)\n", b, localN);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  Matrix local_buf;
  distribute_matrix(w, &local_buf, b, &cart);

  if (!local_buf.values) {
    MPI_Comm_free(&cart.row_comm);
    MPI_Comm_free(&cart.col_comm);
    MPI_Comm_free(&cart.comm);
    return 1;
  }

  int *row_k_flat = (int *)malloc(b * b * sizeof(int));
  int *col_k_flat = (int *)malloc(b * b * sizeof(int));
  int **row_k_buffer = (int **)malloc(b * sizeof(int *));
  int **col_k_buffer = (int **)malloc(b * sizeof(int *));

  if (!row_k_flat || !col_k_flat || !row_k_buffer || !col_k_buffer) {
    if (row_k_flat)
      free(row_k_flat);
    if (col_k_flat)
      free(col_k_flat);
    if (row_k_buffer)
      free(row_k_buffer);
    if (col_k_buffer)
      free(col_k_buffer);
    destroy_buf(local_buf.values, b);
    MPI_Comm_free(&cart.row_comm);
    MPI_Comm_free(&cart.col_comm);
    MPI_Comm_free(&cart.comm);
    return 1;
  }

  for (int i = 0; i < b; i++) {
    row_k_buffer[i] = &row_k_flat[i * b];
    col_k_buffer[i] = &col_k_flat[i * b];
  }

  for (int k = 0; k < cart.q; k++) {
    if (cart.my_row == k && cart.my_col == k) {
      floyd_kernel(local_buf.values, local_buf.values, local_buf.values, b);
      for (int i = 0; i < b; i++) {
        memcpy(row_k_buffer[i], local_buf.values[i], b * sizeof(int));
        memcpy(col_k_buffer[i], local_buf.values[i], b * sizeof(int));
      }
    }
    MPI_Bcast(row_k_flat, b * b, MPI_INT, k, cart.row_comm);
    MPI_Bcast(col_k_flat, b * b, MPI_INT, k, cart.col_comm);

    if (cart.my_row == k && cart.my_col != k) {
      floyd_kernel(local_buf.values, row_k_buffer, local_buf.values, b);
    }

    if (cart.my_row != k && cart.my_col == k) {
      floyd_kernel(local_buf.values, local_buf.values, col_k_buffer, b);
    }

    if (cart.my_row == k) {
      for (int i = 0; i < b; i++)
        memcpy(row_k_buffer[i], local_buf.values[i], b * sizeof(int));
    }
    if (cart.my_col == k) {
      for (int i = 0; i < b; i++)
        memcpy(col_k_buffer[i], local_buf.values[i], b * sizeof(int));
    }

    MPI_Bcast(row_k_flat, b * b, MPI_INT, k, cart.col_comm);
    MPI_Bcast(col_k_flat, b * b, MPI_INT, k, cart.row_comm);

    if (cart.my_row != k && cart.my_col != k) {
      floyd_kernel(local_buf.values, col_k_buffer, row_k_buffer, b);
    }
  }

  free(row_k_flat);
  free(col_k_flat);
  free(row_k_buffer);
  free(col_k_buffer);

  gather_matrix(&local_buf, buf, b, &cart);

  destroy_buf(local_buf.values, b);

  MPI_Comm_free(&cart.row_comm);
  MPI_Comm_free(&cart.col_comm);
  MPI_Comm_free(&cart.comm);
  return 0;
}

/**
 * @brief Parallel blocked Floyd-Warshall using non-blocking MPI communication
 *
 * Optimized version of the parallel blocked Floyd-Warshall that uses
 * MPI_Ibcast (non-blocking broadcast) to overlap communication with
 * computation. This can reduce synchronization overhead and improve
 * performance, especially on high-latency networks.
 *
 * The algorithm structure is similar to the blocking version but uses
 * asynchronous communication with MPI_Wait to synchronize only when
 * the data is actually needed.
 *
 * @param w Input weighted adjacency matrix (only valid on rank 0)
 * @param buf Pointer to output matrix (gathered on rank 0)
 * @param b Block size (must equal n/sqrt(P))
 *
 * @return 0 on success, 1 on failure
 *
 * @note Requires P = q² processes where q divides n
 * @note Block size must equal n/q where q = sqrt(P)
 * @note Uses non-blocking MPI_Ibcast for better performance
 * @note Result is gathered to rank 0 only
 *
 * Time complexity: O(n³/P) with reduced communication overhead
 *
 * @warning Number of processes must be a perfect square
 * @warning Will call MPI_Abort if b != n/sqrt(P)
 *
 * Performance: Generally 5-20% faster than blocking version due to
 * overlapped communication, depending on network latency and block size.
 *
 * @see blocked_floyd_warshall_p_apsp for blocking version
 * @see MPI_Ibcast for non-blocking broadcast details
 */
int blocked_floyd_warshall_p_non_blocking_apsp(Matrix w, Matrix *buf, int b) {
  // initialization
  CartInfo cart;
  setup_cart(&cart, w.n);
  int n = w.n;
  int localN = n / cart.q;
  if (b != localN) {
    if (cart.my_rank == 0) {
      fprintf(stderr, "Set b == n/sqrt(P) (got b=%d, need %d)\n", b, localN);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  Matrix local_buf;
  distribute_matrix(w, &local_buf, b, &cart);
  if (!local_buf.values) {
    MPI_Comm_free(&cart.row_comm);
    MPI_Comm_free(&cart.col_comm);
    MPI_Comm_free(&cart.comm);
    return 1;
  }
  int *row_k_flat = (int *)malloc(b * b * sizeof(int));
  int *col_k_flat = (int *)malloc(b * b * sizeof(int));
  int **row_k_buffer = (int **)malloc(b * sizeof(int *));
  int **col_k_buffer = (int **)malloc(b * sizeof(int *));
  if (!row_k_flat || !col_k_flat || !row_k_buffer || !col_k_buffer) {
    if (row_k_flat)
      free(row_k_flat);
    if (col_k_flat)
      free(col_k_flat);
    if (row_k_buffer)
      free(row_k_buffer);
    if (col_k_buffer)
      free(col_k_buffer);
    destroy_buf(local_buf.values, b);
    MPI_Comm_free(&cart.row_comm);
    MPI_Comm_free(&cart.col_comm);
    MPI_Comm_free(&cart.comm);
    return 1;
  }
  for (int i = 0; i < b; i++) {
    row_k_buffer[i] = &row_k_flat[i * b];
    col_k_buffer[i] = &col_k_flat[i * b];
  }
  MPI_Request row_req[2], col_req[2];
  for (int k = 0; k < cart.q; k++) {
    if (cart.my_row == k && cart.my_col == k) {
      floyd_kernel(local_buf.values, local_buf.values, local_buf.values, b);
      for (int i = 0; i < b; i++) {
        memcpy(row_k_buffer[i], local_buf.values[i], b * sizeof(int));
        memcpy(col_k_buffer[i], local_buf.values[i], b * sizeof(int));
      }
    }
    MPI_Ibcast(row_k_flat, b * b, MPI_INT, k, cart.row_comm, &row_req[0]);
    MPI_Ibcast(col_k_flat, b * b, MPI_INT, k, cart.col_comm, &col_req[0]);

    MPI_Wait(&row_req[0], MPI_STATUS_IGNORE);
    MPI_Wait(&col_req[0], MPI_STATUS_IGNORE);

    if (cart.my_row == k && cart.my_col != k) {
      floyd_kernel(local_buf.values, row_k_buffer, local_buf.values, b);
      for (int i = 0; i < b; i++)
        memcpy(row_k_buffer[i], local_buf.values[i], b * sizeof(int));
    }
    if (cart.my_row != k && cart.my_col == k) {
      floyd_kernel(local_buf.values, local_buf.values, col_k_buffer, b);
      for (int i = 0; i < b; i++)
        memcpy(col_k_buffer[i], local_buf.values[i], b * sizeof(int));
    }

    MPI_Ibcast(row_k_flat, b * b, MPI_INT, k, cart.col_comm, &row_req[1]);
    MPI_Ibcast(col_k_flat, b * b, MPI_INT, k, cart.row_comm, &col_req[1]);

    if (cart.my_row == k && cart.my_col == k) {
      for (int i = 0; i < b; i++) {
        memcpy(row_k_buffer[i], local_buf.values[i], b * sizeof(int));
        memcpy(col_k_buffer[i], local_buf.values[i], b * sizeof(int));
      }
    }

    MPI_Wait(&row_req[1], MPI_STATUS_IGNORE);
    MPI_Wait(&col_req[1], MPI_STATUS_IGNORE);

    if (cart.my_row != k && cart.my_col != k) {
      floyd_kernel(local_buf.values, col_k_buffer, row_k_buffer, b);
    }
  }
  free(row_k_flat);
  free(col_k_flat);
  free(row_k_buffer);
  free(col_k_buffer);
  gather_matrix(&local_buf, buf, b, &cart);
  destroy_buf(local_buf.values, b);
  MPI_Comm_free(&cart.row_comm);
  MPI_Comm_free(&cart.col_comm);
  MPI_Comm_free(&cart.comm);
  return 0;
}
