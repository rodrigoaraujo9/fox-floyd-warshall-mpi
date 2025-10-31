#ifndef TYPES_H
#define TYPES_H

#include <mpi.h>

/**
 * @brief Safe minimum macro for integer comparison.
 *
 * Returns the smaller of two integer values. This macro evaluates each
 * parameter exactly once, making it safe for use with function calls
 * or expressions with side effects.
 *
 * @param a First integer value
 * @param b Second integer value
 * @return The minimum of a and b
 *
 * Example usage:
 * @code
 * int x = 5, y = 3;
 * int min_val = MIN(x, y); // Returns 3
 * @endcode
 */
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/**
 * @brief Matrix structure for n×n integer matrices.
 *
 * Represents a square matrix of integers with dynamic memory allocation.
 * The matrix uses a 2D array representation where values is an array of
 * pointers to row arrays, enabling efficient row-wise access patterns.
 *
 * Memory layout:
 * - values: Array of n pointers to integer arrays
 * - Each values[i] points to an array of n integers representing row i
 * - n: Dimension of the square matrix (n × n)
 *
 * @note Matrices should be allocated using allocate_matrix() and destroyed
 *       with destroy_matrix() to ensure proper memory management.
 * @see allocate_matrix
 * @see destroy_matrix
 */
typedef struct {
    int **values; /**< 2D array of matrix elements, row-major order */
    int n; /**< Dimension of the square matrix (n × n) */
} Matrix;

/**
 * @brief Input/Output file descriptor structure.
 *
 * Simple container for file paths used by I/O functions. This structure
 * holds pointers to input and output file names for reading graph data
 * and writing/verifying results.
 *
 * @note The strings are not copied internally - the structure only holds
 *       pointers to existing string literals or allocated memory.
 * @warning Users must ensure the pointed-to strings remain valid for the
 *          lifetime of the IO_Files structure.
 */
typedef struct {
    char *i; /**< Input file path (graph adjacency matrix) */
    char *o; /**< Output file path (expected shortest path matrix) */
} IO_Files;

/**
 * @brief MPI Cartesian grid information for parallel algorithms.
 *
 * Contains all necessary information for processes operating in a 2D grid
 * topology. This structure is used by parallel blocked Floyd-Warshall
 * implementations to manage process coordination and communication patterns.
 *
 * The grid organizes processes in a q×q arrangement where q = √p, enabling
 * efficient matrix distribution and communication for blocked algorithms.
 *
 * @note The Cartesian communicator and sub-communicators should be freed
 *       with MPI_Comm_free() when no longer needed.
 * @see blocked_floyd_warshall_p_apsp
 * @see blocked_floyd_warshall_p_non_blocking_apsp
 */
typedef struct {
    int p; /**< Total number of processes in MPI_COMM_WORLD */
    MPI_Comm comm; /**< 2D Cartesian communicator for the process grid */
    MPI_Comm row_comm; /**< Sub-communicator for processes in the same row */
    MPI_Comm col_comm; /**< Sub-communicator for processes in the same column */
    int q; /**< Grid dimension (q = √p, number of rows/columns) */
    int my_row; /**< Row coordinate of this process in the grid */
    int my_col; /**< Column coordinate of this process in the grid */
    int my_rank; /**< Rank of this process in MPI_COMM_WORLD */
} CartInfo;

#endif // TYPES_H
