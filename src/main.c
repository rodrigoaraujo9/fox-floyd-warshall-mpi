#include "../includes/algorithms.h"
#include "../includes/blocked_fw.h"
#include "../includes/comm.h"
#include "../includes/io.h"
#include "../includes/matrix.h"
#include "../includes/types.h"
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
  ALG_SLOW,
  ALG_REPEATED_SQUARING,
  ALG_FLOYD_WARSHALL,
  ALG_BLOCKED_FW,
  ALG_BLOCKED_FW_MPI
} Algorithm;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int world_size = 0, world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (argc != 4) {
    if (world_rank == 0) {
      fprintf(stderr,
              "Usage: %s <algorithm> <input_matrix> <expected_output> -> "
              "algorithm {slow, rs , fw, blocked}\n",
              argv[0]);
    }
    MPI_Finalize();
    return 1;
  }

  char *alg_name = argv[1];
  char *input_file = argv[2];
  char *output_file = argv[3];

  // Parse algorithm
  Algorithm alg;
  if (strcmp(alg_name, "slow") == 0) {
    alg = ALG_SLOW;
  } else if (strcmp(alg_name, "rs") == 0) {
    alg = ALG_REPEATED_SQUARING;
  } else if (strcmp(alg_name, "fw") == 0) {
    alg = ALG_FLOYD_WARSHALL;
  } else if (strcmp(alg_name, "blocked") == 0) {
    alg = ALG_BLOCKED_FW;
  } else if (strcmp(alg_name, "mpi") == 0) {
    alg = ALG_BLOCKED_FW_MPI;
  } else {
    if (world_rank == 0) {
      fprintf(stderr, "alg not supported!");
    }
    MPI_Finalize();
    return 1;
  }

  // Read input matrix (all ranks for MPI, rank 0 only for others)
  Matrix input = {NULL, 0};
  if (alg == ALG_BLOCKED_FW_MPI || world_rank == 0) {
    input = read_input_matrix_from_file(input_file);
    if (input.values == NULL) {
      if (world_rank == 0)
        fprintf(stderr, "Failed to read input matrix: %s\n", input_file);
      MPI_Finalize();
      return 1;
    }
  }

  // Read expected output (rank 0 only)
  Matrix expected = {NULL, 0};
  if (world_rank == 0) {
    int n = input.n;
    expected = read_output_matrix_from_file(output_file, n);
    if (expected.values == NULL) {
      fprintf(stderr, "Failed to read expected matrix: %s\n", output_file);
      destroy_matrix(input);
      MPI_Finalize();
      return 1;
    }
  }

  // Run algorithm
  Matrix result = {NULL, 0};
  double t0, t1;
  int rc = 0;

  if (alg == ALG_BLOCKED_FW_MPI) {
    // parallel alg
    CartInfo cart;
    setup_cart(&cart, input.n);
    int b = (int)(input.n / sqrt(cart.p));

    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    rc = blocked_floyd_warshall_p_apsp(input, &result, b);
    t1 = MPI_Wtime();

    if (world_rank == 0) {
      printf("%s,%d,%.6f,%d,%d\n", alg_name, input.n, t1 - t0, world_size, b);
    }
  } else {
    // seq algorythm
    if (world_rank == 0) {
      t0 = MPI_Wtime();

      switch (alg) {
      case ALG_SLOW:
        rc = slow_apsp(input, &result);
        break;
      case ALG_REPEATED_SQUARING:
        rc = repeated_squaring_apsp(input, &result);
        break;
      case ALG_FLOYD_WARSHALL:
        rc = floyd_warshall_apsp(input, &result);
        break;
      case ALG_BLOCKED_FW: {
        int b = (int)sqrt(input.n);
        while (input.n % b != 0)
          b--;
        rc = blocked_floyd_warshall_apsp(input, &result, b);
        break;
      }
      default:
        rc = 1;
      }

      t1 = MPI_Wtime();

      if (rc != 0) {
        fprintf(stderr, "Algorithm failed\n");
      } else {
        printf("%s,%d,%.6f\n", alg_name, input.n, t1 - t0);
      }
    }
  }

  // verify res
  if (world_rank == 0 && rc == 0) {
    assert_apsp(result, expected, alg_name);
  }

  // cleanup
  destroy_matrix(input);
  if (world_rank == 0) {
    destroy_matrix(expected);
  }
  destroy_matrix(result);

  MPI_Finalize();
  return rc;
}
