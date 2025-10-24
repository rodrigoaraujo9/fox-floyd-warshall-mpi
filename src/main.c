#include "../includes/algorithms.h"
#include "../includes/blocked_fw.h"
#include "../includes/io.h"
#include "../includes/matrix.h"
#include "../includes/types.h"
#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

IO_Files io_files[] = {
    {.i = "matrix_examples/input5", .o = "matrix_examples/output5"},
    {.i = "matrix_examples/input6", .o = "matrix_examples/output6"},
    {.i = "matrix_examples/input300", .o = "matrix_examples/output300"},
    {.i = "matrix_examples/input600", .o = "matrix_examples/output600"},
    {.i = "matrix_examples/input900", .o = "matrix_examples/output900"},
    {.i = "matrix_examples/input1200", .o = "matrix_examples/output1200"},
};

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int world_size = 0, world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (argc != 3) {
    if (world_rank == 0) {
      fprintf(stderr, "Usage: %s <input_matrix> <expected_output>\n", argv[0]);
    }
    MPI_Finalize();
    return 1;
  }

  Matrix input = read_input_matrix_from_file(argv[1]);
  if (input.values == NULL) {
    if (world_rank == 0)
      fprintf(stderr, "Failed to read input matrix: %s\n", argv[1]);
    MPI_Finalize();
    return 1;
  }
  Matrix output_file_matrix = read_output_matrix_from_file(argv[2], input.n);
  if (output_file_matrix.values == NULL) {
    if (world_rank == 0)
      fprintf(stderr, "Failed to read expected matrix: %s\n", argv[2]);
    destroy_matrix(input);
    MPI_Finalize();
    return 1;
  }

  CartInfo cart;
  setup_cart(&cart, input.n);
  int p = cart.grid_dim;
  int b = input.n / p;

  Matrix rs_out = (Matrix){0};
  Matrix fw_out = (Matrix){0};

  double t0, t1;

  if (world_rank == 0) {
    t0 = MPI_Wtime();
    if (repeated_squaring_apsp(input, &rs_out) != 0) {
      fprintf(stderr, "Repeated Squaring APSP failed\n");
      destroy_matrix(input);
      destroy_matrix(output_file_matrix);
      MPI_Finalize();
      return 1;
    }
    t1 = MPI_Wtime();
    printf("Repeated Squaring APSP: Time = %.6f s\n", t1 - t0);

    t0 = MPI_Wtime();
    if (floyd_warshall_apsp(input, &fw_out) != 0) {
      fprintf(stderr, "Floyd–Warshall APSP failed\n");
      destroy_matrix(input);
      destroy_matrix(output_file_matrix);
      destroy_matrix(rs_out);
      MPI_Finalize();
      return 1;
    }
    t1 = MPI_Wtime();
    printf("Floyd–Warshall APSP:   Time = %.6f s\n", t1 - t0);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  Matrix mpi_out = (Matrix){0};
  t0 = MPI_Wtime();
  int rc = blocked_floyd_warshall_p_apsp(input, &mpi_out, b);
  t1 = MPI_Wtime();
  if (rc != 0) {
    if (world_rank == 0)
      fprintf(stderr, "MPI Blocked FW APSP failed\n");
    destroy_matrix(input);
    destroy_matrix(output_file_matrix);
    if (world_rank == 0) {
      destroy_matrix(rs_out);
      destroy_matrix(fw_out);
    }
    MPI_Finalize();
    return 1;
  }

  if (world_rank == 0) {
    printf("MPI Blocked FW APSP:   Time = %.6f s (P=%d, p=%d, b=%d)\n", t1 - t0,
           world_size, p, b);

    assert_apsp(mpi_out, output_file_matrix, "MPI Blocked FW APSP");
    assert_apsp(rs_out, output_file_matrix, "Repeated Squaring APSP");
    assert_apsp(fw_out, output_file_matrix, "Floyd–Warshall APSP");
  }

  destroy_matrix(input);
  destroy_matrix(output_file_matrix);
  destroy_matrix(mpi_out);
  if (world_rank == 0) {
    destroy_matrix(rs_out);
    destroy_matrix(fw_out);
  }

  MPI_Finalize();
  return 0;
}
