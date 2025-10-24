#ifndef TYPES_H
#define TYPES_H

#include <mpi.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

typedef struct {
  int **values;
  int n;
} Matrix;

typedef struct {
  char *i;
  char *o;
} IO_Files;

typedef struct {
  int p;             // total number of processes
  int tag;
  MPI_Comm comm;     // 2D Cartesian communicator
  MPI_Comm row_comm; // subcomm for my row (broadcast A in Fox)
  MPI_Comm col_comm; // subcomm for my column (optional sync/shift of B)
  int q;             // sqrt(p)
  int my_row;        // row coordinate in the grid
  int my_col;        // col coordinate in the grid
  int my_rank;       // rank in MPI_COMM_WORLD
} CartInfo;

#endif // TYPES_H
