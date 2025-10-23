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
  MPI_Comm cart_comm;
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  int p_row;
  int p_col;
  int g_rows;
  int g_cols;
  int my_rank;
  int cart_rank;
  int row_rank;
  int size;
} CartInfo;

#endif // TYPES_H
