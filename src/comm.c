#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Sets up 2D Cartesian process topology for matrix operations
 *
 * Creates a 2D Cartesian grid of processes (sqrt(P) × sqrt(P)) and
 * establishes row and column communicators for efficient collective
 * operations. Each process is assigned coordinates (row, col) in the grid.
 *
 * The function performs validation to ensure:
 * - Number of processes is a perfect square
 * - Matrix size is divisible by sqrt(P)
 *
 * @param cart Pointer to CartInfo structure to be initialized
 * @param n Matrix dimension (must be divisible by sqrt(P))
 *
 * @note Calls MPI_Abort on failure (entire job terminates)
 * @note Creates three communicators that must be freed by caller:
 *       - cart->comm (2D Cartesian)
 *       - cart->row_comm (row subcommunicator)
 *       - cart->col_comm (column subcommunicator)
 *
 * @warning Caller must free communicators with MPI_Comm_free
 * @warning Aborts if P is not a perfect square or n % sqrt(P) != 0
 *
 * CartInfo fields populated:
 * - p: Total number of processes
 * - q: Grid dimension (sqrt(P))
 * - my_rank: Process rank in MPI_COMM_WORLD
 * - my_row: Process row coordinate (0 to q-1)
 * - my_col: Process column coordinate (0 to q-1)
 * - comm: 2D Cartesian communicator
 * - row_comm: Communicator for processes in same row
 * - col_comm: Communicator for processes in same column
 *
 * Example topology for P=4:
 * @code
 * Process grid:     Row comms:      Col comms:
 * (0,0) (0,1)       {0,1}           {0,2}
 * (1,0) (1,1)       {2,3}           {1,3}
 * @endcode
 *
 * @see blocked_floyd_warshall_p_apsp for usage example
 * @see MPI_Cart_create, MPI_Cart_sub
 */
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

/**
 * @brief Distributes a matrix from root process to all processes in blocks
 *
 * Scatters the global matrix w from rank 0 to all processes in the 2D grid.
 * Each process receives a b × b block corresponding to its position in the
 * grid. Process at (row, col) receives the block starting at global position
 * (row*b, col*b).
 *
 * Communication pattern:
 * - Rank 0: Copies its own block directly and sends blocks to all others
 * - Other ranks: Receive their respective blocks from rank 0
 *
 * @param w Global input matrix (only valid on rank 0, unused on others)
 * @param local_buf Pointer to Matrix structure for storing local block
 * @param b Block size (dimension of each local block)
 * @param cart Pointer to CartInfo with topology information
 *
 * @note local_buf->values is allocated by this function
 * @note local_buf->n is set to the global matrix size n (not block size b)
 * @note Uses point-to-point communication (MPI_Send/MPI_Recv)
 *
 * @warning Caller must free local_buf->values after use
 * @warning w must be allocated and valid on rank 0
 *
 * Memory layout:
 * - Global matrix w: n × n on rank 0
 * - Local blocks: b × b on each rank where b = n/q
 *
 * Example for 4×4 matrix distributed to 4 processes (b=2):
 * @code
 * Global matrix w:     Process blocks:
 * [A B | C D]          Rank 0: [A B]   Rank 1: [C D]
 * [E F | G H]                  [E F]           [G H]
 * ----+----
 * [I J | K L]          Rank 2: [I J]   Rank 3: [K L]
 * [M N | O P]                  [M N]           [O P]
 * @endcode
 *
 * @see gather_matrix for the inverse operation
 * @see setup_cart for CartInfo initialization
 */
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

/**
 * @brief Gathers matrix blocks from all processes to root process
 *
 * Collects b × b blocks from all processes and assembles them into a
 * complete n × n matrix on rank 0. This is the inverse operation of
 * distribute_matrix.
 *
 * Communication pattern:
 * - Rank 0: Receives blocks from all other processes and copies its own
 * - Other ranks: Send their local blocks to rank 0
 *
 * @param local_buf Pointer to local block matrix (b × b on each process)
 * @param global_buf Pointer to global matrix (allocated only on rank 0)
 * @param b Block size (dimension of each local block)
 * @param cart Pointer to CartInfo with topology information
 *
 * @note global_buf->values is allocated only on rank 0
 * @note global_buf is uninitialized on non-root processes
 * @note local_buf->n contains global matrix size n
 *
 * @warning Caller must free global_buf->values on rank 0
 * @warning global_buf is only valid on rank 0 after this call
 *
 * Assembly process:
 * Blocks are gathered in row-major order of the process grid and placed
 * at their corresponding positions in the global matrix.
 *
 * Example for 4 processes (b=2) assembling 4×4 matrix:
 * @code
 * Process blocks:      Assembled global matrix:
 * Rank 0: [A B]        [A B | C D]
 *         [E F]        [E F | G H]
 * Rank 1: [C D]        ----+----
 *         [G H]        [I J | K L]
 * Rank 2: [I J]        [M N | O P]
 *         [M N]
 * Rank 3: [K L]
 *         [O P]
 * @endcode
 *
 * Time complexity: O(n²/P) send operations per process, sequential assembly
 *
 * @see distribute_matrix for the inverse operation
 * @see setup_cart for CartInfo initialization
 */
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
