#include "../includes/blocked_fw.h"
#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void setup_cart(CartInfo *cart, int n) {
  MPI_Comm_size(MPI_COMM_WORLD, &cart->p);
  MPI_Comm_rank(MPI_COMM_WORLD, &cart->my_rank);

  cart->q = (int)(sqrt((double)cart->p));
  if (cart->q * cart->q != cart->p) {
    if (cart->my_rank == 0)
      fprintf(stderr, "P must be a perfect square (got %d).\n", cart->p);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (n % cart->q != 0) {
    if (cart->my_rank == 0)
      fprintf(stderr, "Matrix size n=%d must be divisible by sqrt(P)=%d.\n", n,
              cart->q);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int dims[2] = {cart->q, cart->q};
  int periods[2] = {0, 0};
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart->comm);
  if (cart->comm == MPI_COMM_NULL) {
    if (cart->my_rank == 0)
      fprintf(stderr, "MPI_Cart_create failed.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int comm_rank;
  MPI_Comm_rank(cart->comm, &comm_rank);
  int coords[2];
  MPI_Cart_coords(cart->comm, comm_rank, 2, coords);
  cart->my_row = coords[0];
  cart->my_col = coords[1];

  int remain[2] = {0, 1};
  MPI_Cart_sub(cart->comm, remain, &cart->row_comm);

  remain[0] = 1;
  remain[1] = 0;
  MPI_Cart_sub(cart->comm, remain, &cart->col_comm);
}

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

int broadcast_block_to_row(int **a, int **buf, int na, int ma, int step,
                           CartInfo *cart) {
  int root;
  int count;

  count = ma * na / cart->p;

  if (cart->my_rank ==
      cart->my_row * cart->q + (cart->my_row + step) % cart->q) {
    memcpy(buf, a, count * sizeof(float));
  }

  root = (cart->my_row + step % cart->q) % cart->q;
  MPI_Bcast(buf, count, MPI_INT, root, cart->row_comm);

  return 0;
}

int circular_shift(int **b, int mb, int nb, CartInfo *cart) {
  int dest;
  int source;
  int tag = 0;
  MPI_Status status;

  source = (cart->q + cart->my_row + 1) % cart->q;
  dest = (cart->q + cart->my_row - 1) % cart->q;

  MPI_Sendrecv_replace(b, mb * nb / cart->p, MPI_FLOAT, dest, tag, source, tag,
                       cart->col_comm, &status);

  return 0;
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
  int localN = n / cart.p;

  if (b != localN) {
    if (cart.my_rank == 0) {
      fprintf(stderr, "Set b == n/sqrt(P) (got b=%d, need %d)\n", b, localN);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  copy_matrix(w, buf);
  if (!buf->values) {
    MPI_Comm_free(&cart.row_comm);
    MPI_Comm_free(&cart.col_comm);
    MPI_Comm_free(&cart.comm);
    return 1;
  }

  int *wkk_buf = (int *)malloc(b * b * sizeof(int));
  int *wkj_buf = (int *)malloc(b * b * sizeof(int));
  int *wik_buf = (int *)malloc(b * b * sizeof(int));

  int **wkk_ptr = (int **)malloc(b * sizeof(int *));
  int **wkj_ptr = (int **)malloc(b * sizeof(int *));
  int **wik_ptr = (int **)malloc(b * sizeof(int *));
  for (int i = 0; i < b; i++) {
    wkk_ptr[i] = &wkk_buf[i * b];
    wkj_ptr[i] = &wkj_buf[i * b];
    wik_ptr[i] = &wik_buf[i * b];
  }

  int B = n / b;

  // computation
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

  free(wkk_buf);
  free(wkj_buf);
  free(wik_buf);
  free(wkk_ptr);
  free(wkj_ptr);
  free(wik_ptr);

  MPI_Comm_free(&cart.row_comm);
  MPI_Comm_free(&cart.col_comm);
  MPI_Comm_free(&cart.comm);
  return 0;
}
