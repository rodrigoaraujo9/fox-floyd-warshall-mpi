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

void distribute_matrix(Matrix w, Matrix *local_buf, int b, CartInfo *cart) {
  int n = w.n;

  // metadata - each block needs to know where it fits in the global cart
  local_buf->n = n;

  // though it's size is the actual block size
  local_buf->values = allocate_matrix(b);

  if (cart->my_rank == 0) {
    // distribute the matrix
    for (int row = 0; row < cart->q; row++) {
      for (int col = 0; col < cart->q; col++) {
        int dst_rank = row * cart->q + col;
        for (int i = 0; i < b; i++) {
          int glob_row = row * b + i;
          int col_init = col * b;
          if (dst_rank == 0) {
            memcpy(local_buf->values[i], &w.values[glob_row][col_init],
                   b * sizeof(int));
          } else {
            MPI_Send(&w.values[glob_row][col_init], b, MPI_INT, dst_rank, i,
                     cart->comm);
          }
        }
      }
    }
  } else {
    // receive local block
    for (int i = 0; i < b; i++) {
      MPI_Recv(local_buf->values[i], b, MPI_INT, 0, i, cart->comm,
               MPI_STATUS_IGNORE);
    }
  }
}

void gather_matrix(Matrix *local_buf, Matrix *global_buf, int b,
                   CartInfo *cart) {
  int n = local_buf->n;
  if (cart->my_rank == 0) {
    // gather each local block and assemble result matrix
    global_buf->n = n;
    global_buf->values = allocate_matrix(n);

    // iterate over process [row][col]
    for (int row = 0; row < cart->q; row++) {
      for (int col = 0; col < cart->q; col++) {

        // rank in global communicator
        int src_rank = row * cart->q + col;

        for (int i = 0; i < b; i++) {
          int glob_row = row * b + i;
          int col_init = col * b;
          if (src_rank == 0) {
            memcpy(&global_buf->values[glob_row][col_init],
                   local_buf->values[i], b * sizeof(int));
          } else {
            MPI_Recv(&global_buf->values[glob_row][col_init], b, MPI_INT,
                     src_rank, i, cart->comm, MPI_STATUS_IGNORE);
          }
        }
      }
    }

  } else {
    // send local block to process with rank 0
    for (int i = 0; i < b; i++) {
      MPI_Send(local_buf->values[i], b, MPI_INT, 0, i, cart->comm);
    }
  }
}
