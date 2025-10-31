#include "../includes/algorithms.h"
#include "../includes/blocked_fw.h"
#include "../includes/io.h"
#include "../includes/matrix.h"
#include "../includes/types.h"
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

/**
 * @enum Algorithm
 * @brief Enumeration of supported APSP algorithm variants.
 */
typedef enum {
    ALG_SLOW, /**< Naive O(n⁴) algorithm */
    ALG_REPEATED_SQUARING, /**< Repeated squaring (matrix multiplication–based) */
    ALG_FLOYD_WARSHALL, /**< Classic Floyd–Warshall algorithm */
    ALG_BLOCKED_FW, /**< Cache-blocked Floyd–Warshall algorithm */
    ALG_BLOCKED_FW_MPI, /**< Parallel blocked Floyd–Warshall using MPI (blocking) */
    ALG_BLOCKED_FW_MPI_NON_BLOCKING /**< Parallel blocked Floyd–Warshall using MPI (non-blocking) */
} Algorithm;

/**
 * @brief Creates the output directory for result logs if it doesn't exist.
 *
 * Ensures that the results/ directory is present before attempting to write CSV
 * benchmark logs. If the directory does not exist, it is created with permissions 0755.
 *
 * @note Does nothing if the directory already exists.
 * @warning Fails silently if directory creation is not permitted.
 * @see log_result_to_file
 */
void create_results_dir() {
  struct stat st = {0};
  if (stat("results", &st) == -1) {
    mkdir("results", 0755);
  }
}

/**
 * @brief Appends a benchmark result entry to a CSV file.
 *
 * Records algorithm execution results including timestamp, algorithm name,
 * matrix size, elapsed time, number of MPI processes, and block size.
 *
 * @param alg_name Algorithm name (e.g., "fw", "mpi", etc.)
 * @param n Matrix size (dimension n × n)
 * @param exec_time Execution time in seconds
 * @param world_size Number of MPI processes (1 for sequential)
 * @param block_size Block size used (0 if not applicable)
 *
 * @note Automatically creates results/ directory if missing.
 * @note If the file is new, a CSV header is written.
 * @warning File I/O errors are ignored silently.
 */
void log_result_to_file(const char *alg_name, int n, double exec_time,
                        int world_size, int block_size) {
  char filename[256];
  time_t now;
  time(&now);
  struct tm *t = localtime(&now);

  create_results_dir();

  if (world_size > 1) {
    snprintf(filename, sizeof(filename), "results/%s_n%d_p%d_b%d.csv", alg_name,
             n, world_size, block_size);
  } else {
    snprintf(filename, sizeof(filename), "results/%s_n%d.csv", alg_name, n);
  }

  FILE *fp = fopen(filename, "a");
  if (fp) {
    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);

    if (size == 0) {
      // Write header
      if (world_size > 1) {
        fprintf(
            fp,
            "timestamp,algorithm,matrix_size,time_sec,processes,block_size\n");
      } else {
        fprintf(fp, "timestamp,algorithm,matrix_size,time_sec\n");
      }
    }

    char timestamp[64];
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", t);

    if (world_size > 1) {
      fprintf(fp, "%s,%s,%d,%.6f,%d,%d\n", timestamp, alg_name, n, exec_time,
              world_size, block_size);
    } else {
      fprintf(fp, "%s,%s,%d,%.6f\n", timestamp, alg_name, n, exec_time);
    }

    fclose(fp);
  }
}

/**
 * @brief Prints a formatted summary of benchmark results to stdout.
 *
 * Displays algorithm configuration, execution time, performance metrics
 * (GFLOPs), and verification results in a readable ASCII banner.
 *
 * @param alg_name Algorithm name
 * @param n Matrix dimension
 * @param exec_time Execution time in seconds
 * @param world_size Number of MPI processes
 * @param block_size Block size used (if applicable)
 * @param passed 1 if verification succeeded, 0 otherwise
 *
 * @note Also prints efficiency per process for parallel runs.
 * @see log_result_to_file
 */
void print_result_banner(const char *alg_name, int n, double exec_time,
                         int world_size, int block_size, int passed) {
  printf("\n");
  printf("═══════════════════════════════════════════════════════════\n");
  printf("  Algorithm:     %s\n", alg_name);
  printf("  Matrix Size:   %d × %d\n", n, n);
  printf("  Execution:     %.6f seconds\n", exec_time);

  if (world_size > 1) {
    printf("  Processes:     %d\n", world_size);
    printf("  Block Size:    %d\n", block_size);
  }

  printf("  Verification:  %s\n", passed ? "PASSED" : "FAILED");
  printf("═══════════════════════════════════════════════════════════\n");

  // Performance metrics
  long long operations = (long long)n * n * n;
  double gflops = operations / (exec_time * 1e9);
  printf("  Operations:    %lld\n", operations);
  printf("  Performance:   %.3f GFLOP/s\n", gflops);

  if (world_size > 1) {
    double parallel_efficiency = gflops / world_size;
    printf("  Efficiency:    %.3f GFLOP/s per process\n", parallel_efficiency);
  }

  printf("═══════════════════════════════════════════════════════════\n");

  printf("  Time/Process:  %.6f seconds\n", exec_time);

  printf("═══════════════════════════════════════════════════════════\n");
  printf("\n\n");
  printf("Result logged to: results/%s_n%d", alg_name, n);
  if (world_size > 1) {
    printf("_p%d_b%d", world_size, block_size);
  }
  printf(".csv\n");
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int world_size = 0, world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (argc != 4) {
    if (world_rank == 0) {
      fprintf(stderr, "\n");
      fprintf(stderr,
              "Usage: %s <algorithm> <input_matrix> <expected_output>\n",
              argv[0]);
      fprintf(stderr, "\n");
      fprintf(stderr, "Algorithms:\n");
      fprintf(stderr, "  slow    - Naive O(n^4) implementation\n");
      fprintf(stderr, "  rs      - Repeated squaring\n");
      fprintf(stderr, "  fw      - Floyd-Warshall\n");
      fprintf(stderr, "  blocked - Blocked Floyd-Warshall\n");
      fprintf(stderr, "  mpi     - Parallel blocked FW (MPI)\n");
      fprintf(stderr, "  mpi-nb  - Non-blocking parallel FW (MPI)\n");
      fprintf(stderr, "\n");
    }
    MPI_Finalize();
    return 1;
  }

  char *alg_name = argv[1];
  char *input_file = argv[2];
  char *output_file = argv[3];

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
  } else if (strcmp(alg_name, "mpi-nb") == 0) {
    alg = ALG_BLOCKED_FW_MPI_NON_BLOCKING;
  } else {
    if (world_rank == 0) {
      fprintf(stderr, "Error: Algorithm '%s' not supported!\n", alg_name);
    }
    MPI_Finalize();
    return 1;
  }

  if (world_rank == 0) {
    printf("\n");
    printf("Starting APSP computation...\n");
    printf("Algorithm: %s\n", alg_name);
    if (world_size > 1) {
      printf("Processes: %d\n", world_size);
    }
    printf("\n");
  }

  Matrix input = {NULL, 0};
  if (alg == ALG_BLOCKED_FW_MPI || alg == ALG_BLOCKED_FW_MPI_NON_BLOCKING ||
      world_rank == 0) {
    input = read_input_matrix_from_file(input_file);
    if (input.values == NULL) {
      if (world_rank == 0)
        fprintf(stderr, "Error: Failed to read input matrix: %s\n", input_file);
      MPI_Finalize();
      return 1;
    }
  }

  Matrix expected = {NULL, 0};
  if (world_rank == 0) {
    int n = input.n;
    expected = read_output_matrix_from_file(output_file, n);
    if (expected.values == NULL) {
      fprintf(stderr, "Error: Failed to read expected matrix: %s\n",
              output_file);
      destroy_matrix(input);
      MPI_Finalize();
      return 1;
    }
  }

  Matrix result = {NULL, 0};
  double t0, t1;
  int rc = 0;
  int block_size = 0;

  if (alg == ALG_BLOCKED_FW_MPI) {
    block_size = (int)(input.n / sqrt(world_size));
    flush_cache();
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    rc = blocked_floyd_warshall_p_apsp(input, &result, block_size);
    t1 = MPI_Wtime();
  } else if (alg == ALG_BLOCKED_FW_MPI_NON_BLOCKING) {
    block_size = (int)(input.n / sqrt(world_size));
    flush_cache();
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    rc = blocked_floyd_warshall_p_non_blocking_apsp(input, &result, block_size);
    t1 = MPI_Wtime();
  } else {
    if (world_rank == 0) {
      flush_cache();
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
        block_size = (int)sqrt(input.n);
        while (input.n % block_size != 0)
          block_size--;
        rc = blocked_floyd_warshall_apsp(input, &result, block_size);
        break;
      }
      default:
        rc = 1;
      }
      t1 = MPI_Wtime();

      if (rc != 0) {
        fprintf(stderr, "Error: Algorithm failed\n");
      }
    }
  }

  int passed = 1;
  if (world_rank == 0 && rc == 0) {
    assert_apsp(result, expected, alg_name);
  }

  if (world_rank == 0 && rc == 0) {
    double elapsed = t1 - t0;
    print_result_banner(alg_name, input.n, elapsed, world_size, block_size,
                        passed);
    log_result_to_file(alg_name, input.n, elapsed, world_size, block_size);
  }

  destroy_matrix(input);
  if (world_rank == 0) {
    destroy_matrix(expected);
  }
  destroy_matrix(result);

  MPI_Finalize();
  return rc;
}
