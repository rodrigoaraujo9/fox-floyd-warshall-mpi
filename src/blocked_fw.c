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

  // int B = cart.q;

  int **row_k_buffer = allocate_matrix(b);
  int **col_k_buffer = allocate_matrix(b);

  if (!row_k_buffer || !col_k_buffer) {
    if (row_k_buffer)
      destroy_buf(row_k_buffer, b);
    if (col_k_buffer)
      destroy_buf(col_k_buffer, b);
    destroy_buf(local_buf.values, b);
    MPI_Comm_free(&cart.row_comm);
    MPI_Comm_free(&cart.col_comm);
    MPI_Comm_free(&cart.comm);
    return 1;
  }

  // computing and communication done here

  gather_matrix(&local_buf, buf, b, &cart);

  // currently here to verify functional correctness
  blocked_floyd_warshall_apsp(w, buf, b);

  destroy_buf(row_k_buffer, b);
  destroy_buf(col_k_buffer, b);
  destroy_buf(local_buf.values, b);

  MPI_Comm_free(&cart.row_comm);
  MPI_Comm_free(&cart.col_comm);
  MPI_Comm_free(&cart.comm);
  return 0;
}
