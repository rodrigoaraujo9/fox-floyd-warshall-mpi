#include "../includes/blocked_fw.h"
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

  cart->tag = 1;

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

int send_block_p(Matrix m, int n, CartInfo cart) {
    int b = n * n / cart.p;
    int sqrt_b = (int) sqrt(b);

    int process_rank_count = 1;

    int *send_buf = (int*) malloc(sizeof(int) * b);

    for (int ki = sqrt_b; ki < n; ki += sqrt_b) {
        for (int kj = sqrt_b; kj < n; kj += sqrt_b) {
            int send_buf_count = 0;

            for (int i = ki; i < ki + sqrt_b; i++) {
                for (int j = kj; j < kj + sqrt_b; j++) {
                    send_buf[send_buf_count] = m.values[i][j];
                }
            }

            MPI_Send(send_buf, b * sizeof(int), MPI_INT, process_rank_count, cart.tag, cart.comm);
        }
    }

    free(send_buf);
    return 0;
}

int blocked_floyd_warshall_p_apsp(Matrix w, Matrix *buf, int b) {
    return 0;
}
