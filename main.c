#include <assert.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

IO_Files io_files[] = {
    {.i = "matrix_examples/input5", .o = "matrix_examples/output5"},
    {.i = "matrix_examples/input6", .o = "matrix_examples/output6"},
    {.i = "matrix_examples/input300", .o = "matrix_examples/output300"},
    {.i = "matrix_examples/input600", .o = "matrix_examples/output600"},
    /*
    {.i = "matrix_examples/input900", .o = "matrix_examples/output900"},
    {.i = "matrix_examples/input1200", .o = "matrix_examples/output1200"},
     */
};

Matrix read_input_matrix_from_file(char *file_path);
Matrix read_output_matrix_from_file(char *file_path, int n);

void print_matrix(Matrix matrix, char *title);
void copy_matrix(Matrix matrix_to_copy, Matrix *buf);
void destroy_matrix(Matrix matrix);

int special_matrix_mul(Matrix a, Matrix b, Matrix *buf);

int slow_apsp(Matrix w, Matrix *buf);
int repeated_squaring_apsp(Matrix w, Matrix *buf);
int floyd_warshall_apsp(Matrix w, Matrix *buf);

void assert_apsp(Matrix a, Matrix b, char *title);

void set_block(int **matrix, int block_row, int block_col, int b, int n,
               int **block);

void floyd(int **matrix, int **C, int **A, int **B, int b, int n);

int blocked_floyd_warshall_apsp(Matrix w, Matrix *buf, int b);

void setup_cart(CartInfo *cart, int n);

int blocked_floyd_warshall_p_apsp(Matrix w, Matrix *buf, int b);

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
  Matrix expected = read_output_matrix_from_file(argv[2], input.n);
  if (expected.values == NULL) {
    if (world_rank == 0)
      fprintf(stderr, "Failed to read expected matrix: %s\n", argv[2]);
    destroy_matrix(input);
    MPI_Finalize();
    return 1;
  }

  CartInfo cart;
  setup_cart(&cart, input.n);
  int p = cart.g_rows;
  int b = input.n / p;

  Matrix rs_out = (Matrix){0};
  Matrix fw_out = (Matrix){0};

  double t0, t1;

  if (world_rank == 0) {
    t0 = MPI_Wtime();
    if (repeated_squaring_apsp(input, &rs_out) != 0) {
      fprintf(stderr, "Repeated Squaring APSP failed\n");
      destroy_matrix(input);
      destroy_matrix(expected);
      MPI_Finalize();
      return 1;
    }
    t1 = MPI_Wtime();
    printf("Repeated Squaring APSP: Time = %.6f s\n", t1 - t0);

    t0 = MPI_Wtime();
    if (floyd_warshall_apsp(input, &fw_out) != 0) {
      fprintf(stderr, "Floyd–Warshall APSP failed\n");
      destroy_matrix(input);
      destroy_matrix(expected);
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
    destroy_matrix(expected);
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
  }

  if (world_rank == 0) {
    assert_apsp(mpi_out, expected, "MPI Blocked FW vs Expected");
    assert_apsp(rs_out, expected, "Repeated Squaring vs Expected");
    assert_apsp(fw_out, expected, "Floyd–Warshall vs Expected");
    assert_apsp(mpi_out, fw_out, "MPI Blocked FW vs Floyd–Warshall");
  }

  destroy_matrix(input);
  destroy_matrix(expected);
  destroy_matrix(mpi_out);
  if (world_rank == 0) {
    destroy_matrix(rs_out);
    destroy_matrix(fw_out);
  }

  MPI_Comm_free(&cart.row_comm);
  MPI_Comm_free(&cart.col_comm);
  MPI_Comm_free(&cart.cart_comm);

  MPI_Finalize();
  return 0;
}

void assert_apsp(Matrix a, Matrix b, char *title) {
  assert(a.n == b.n);
  for (int i = 0; i < a.n; i++) {
    for (int j = 0; j < a.n; j++) {
      assert(a.values[i][j] == b.values[i][j]);
    }
  }
  printf("%s: Assertion for matrix of size %d was successful!\n", title, a.n);
}

Matrix read_input_matrix_from_file(char *file_path) {
  FILE *file;
  int n, **values;

  // open the file
  if ((file = fopen(file_path, "r")) == NULL) {
    return (Matrix){values, n};
  }

  // scan the size of the rows and columns
  fscanf(file, "%d", &n);

  // allocate space for each column of the matrix
  values = (int **)malloc(sizeof(int *) * n);

  for (int i = 0; i < n; i++) {
    // allocate space for each line of the matrix
    values[i] = (int *)malloc(sizeof(int) * n);
    for (int j = 0; j < n; j++) {
      // read the value and assign it to the correct position in the matrix
      fscanf(file, "%d", &values[i][j]);
      if (i != j && values[i][j] == 0) {
        values[i][j] = INT_MAX;
      }
    }
  }

  fclose(file);

  return (Matrix){values, n};
}

Matrix read_output_matrix_from_file(char *file_path, int n) {
  FILE *file;
  int **values;

  // open the file
  if ((file = fopen(file_path, "r")) == NULL) {
    return (Matrix){values, n};
  }

  // allocate space for each column of the matrix
  values = (int **)malloc(sizeof(int *) * n);

  for (int i = 0; i < n; i++) {
    // allocate space for each line of the matrix
    values[i] = (int *)malloc(sizeof(int) * n);
    for (int j = 0; j < n; j++) {
      // read the value and assign it to the correct position in the matrix
      fscanf(file, "%d", &values[i][j]);
      if (i != j && values[i][j] == 0) {
        values[i][j] = INT_MAX;
      }
    }
  }

  fclose(file);

  return (Matrix){values, n};
}

void print_matrix(Matrix matrix, char *title) {
  printf("%s (%dx%d)\n", title, matrix.n, matrix.n);
  for (int i = 0; i < matrix.n; i++) {
    for (int j = 0; j < matrix.n; j++) {
      if (matrix.values[i][j] == INT_MAX) {
        printf("0 ");
        continue;
      }
      printf("%d ", matrix.values[i][j]);
    }
    printf("\n");
  }
}

void copy_matrix(Matrix matrix_to_copy, Matrix *buf) {
  (*buf).n = matrix_to_copy.n;
  (*buf).values = (int **)malloc(sizeof(int *) * matrix_to_copy.n);

  for (int i = 0; i < matrix_to_copy.n; i++) {
    (*buf).values[i] = (int *)malloc(sizeof(int) * matrix_to_copy.n);
    for (int j = 0; j < matrix_to_copy.n; j++) {
      (*buf).values[i][j] = matrix_to_copy.values[i][j];
    }
  }
}

void destroy_matrix(Matrix matrix) {
  if (matrix.values != NULL) {
    for (int i = 0; i < matrix.n; i++) {
      if (matrix.values[i] != NULL)
        free(matrix.values[i]);
    }
    free(matrix.values);
  }
}

// given algorythm for matrix mul using dynamic programming
int special_matrix_mul(Matrix a, Matrix b, Matrix *buf) {
  if (a.n != b.n) {
    fprintf(stderr, "matrix sizes must be same: %d != %d\n", a.n, b.n);
    return 1;
  }

  *buf = (Matrix){(int **)malloc(sizeof(int *) * a.n), a.n};
  int help;

  for (int i = 0; i < (*buf).n; i++) {
    (*buf).values[i] = (int *)malloc(sizeof(int) * (*buf).n);
    for (int j = 0; j < (*buf).n; j++) {
      (*buf).values[i][j] = INT_MAX;
      for (int k = 0; k < (*buf).n; k++) {
        if (a.values[i][k] != INT_MAX && b.values[k][j] != INT_MAX) {
          help = (a.values[i][k] + b.values[k][j]);
          if ((*buf).values[i][j] > help) {
            (*buf).values[i][j] = help;
          }
        }
      }
    }
  }

  return 0;
}

// 3 implementation sequential
int slow_apsp(Matrix w, Matrix *buf) {
  Matrix d_next;

  copy_matrix(w, buf);

  // from d(2) to d(n-1)
  for (int k = 2; k < w.n - 1; k++) {
    special_matrix_mul(*buf, w, &d_next);
    if (d_next.values == NULL) {
      fprintf(stderr, "special_matrix_mul failed at k=%d\n", k);
      return 1;
    }
    destroy_matrix(*buf);
    *buf = d_next;
  }
  return 0;
}

// 3.1 implementation sequential
int repeated_squaring_apsp(Matrix w, Matrix *buf) {
  int m;
  Matrix d_2m;
  copy_matrix(w, buf);

  if ((*buf).values == NULL) {
    fprintf(stderr, "Failed to copy initial matrix\n");
    return 1;
  }

  m = 1;

  while (m < w.n - 1) {
    special_matrix_mul((*buf), (*buf), &d_2m);
    if (d_2m.values == NULL) {
      destroy_matrix(*buf);
      fprintf(stderr, "special_matrix_mul failed at m=%d\n", m);
      return 1;
    }

    m *= 2;
    destroy_matrix(*buf);
    (*buf) = d_2m;
  }
  return 0;
}

int floyd_warshall_apsp(Matrix w, Matrix *buf) {
  int new_d, *ri, *rk, rik;
  copy_matrix(w, buf);
  if ((*buf).values == NULL) {
    fprintf(stderr, "Failed to copy initial matrix\n");
    return 1;
  }

  for (int k = 0; k < w.n; k++) {
    rk = (*buf).values[k];
    for (int i = 0; i < w.n; i++) {
      ri = (*buf).values[i];
      rik = ri[k];
      if (rik == INT_MAX)
        continue;
      for (int j = 0; j < w.n; j++) {
        if (rk[j] != INT_MAX) {
          new_d = rik + rk[j];
          if (ri[j] > new_d) {
            ri[j] = new_d;
          }
        }
      }
    }
  }
  return 0;
}

int **get_block(int **matrix, int b_row, int b_col, int b, int n) {
  int r0 = b_row * b;
  int c0 = b_col * b;

  assert(r0 >= 0 && c0 >= 0);
  assert(r0 + b <= n && c0 + b <= n);

  int **block = (int **)malloc(sizeof(int *) * b);
  for (int i = 0; i < b; i++) {
    block[i] = &matrix[r0 + i][c0];
  }
  return block;
}

void floyd(int **matrix, int **C, int **A, int **B, int b, int n) {
  int a_val, b_val, sum;
  for (int k = 0; k < b; k++) {
    for (int i = 0; i < b; i++) {
      a_val = A[i][k];
      if (a_val == INT_MAX)
        continue;
      for (int j = 0; j < b; j++) {
        b_val = B[k][j];
        if (b_val != INT_MAX) {
          sum = a_val + b_val;
          if (C[i][j] > sum) {
            C[i][j] = sum;
          }
        }
      }
    }
  }
}

// for perfect squares only for now
int blocked_floyd_warshall_apsp(Matrix w, Matrix *buf, int b) {
  int n = w.n;
  if (n % b != 0) {
    fprintf(stderr, "Block size must divide matrix size %d\n", n);
    return 1;
  }

  copy_matrix(w, buf);
  if ((*buf).values == NULL) {
    fprintf(stderr, "Failed to copy initial matrix\n");
    return 1;
  }

  int B = n / b;
  int **wkk, **wkj, **wik, **wij;

  for (int k = 0; k < B; k++) {
    // dependant phase
    wkk = get_block((*buf).values, k, k, b, n);
    floyd((*buf).values, wkk, wkk, wkk, b, n);
    // partially dependant phase
    for (int j = 0; j < B; j++) {
      if (j == k)
        continue;
      wkj = get_block((*buf).values, k, j, b, n);
      wkk = get_block((*buf).values, k, k, b, n);
      floyd((*buf).values, wkj, wkk, wkj, b, n);
      free(wkj);
      free(wkk);
    }
    for (int i = 0; i < B; i++) {
      if (i == k)
        continue;
      wik = get_block((*buf).values, i, k, b, n);
      wkk = get_block((*buf).values, k, k, b, n);
      floyd((*buf).values, wik, wik, wkk, b, n);
      free(wkk);

      // independant phase
      for (int j = 0; j < B; j++) {
        if (j == k)
          continue;
        wkj = get_block((*buf).values, k, j, b, n);
        wij = get_block((*buf).values, i, j, b, n);
        floyd((*buf).values, wij, wik, wkj, b, n);
        free(wkj);
        free(wij);
      }
      free(wik);
    }
  }

  return 0;
}

void setup_cart(CartInfo *cart, int n) {
  MPI_Comm_size(MPI_COMM_WORLD, &cart->size);
  MPI_Comm_rank(MPI_COMM_WORLD, &cart->my_rank);

  cart->g_rows = (int)sqrt((double)cart->size);
  cart->g_cols = cart->g_rows;

  if (cart->g_rows * cart->g_cols != cart->size) {
    if (cart->my_rank == 0) {
      fprintf(stderr, "Number of processes must be a perfect square (got %d)\n",
              cart->size);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (n % cart->g_rows != 0) {
    if (cart->my_rank == 0) {
      fprintf(stderr, "Matrix size n=%d must be divisible by sqrt(P)=%d\n", n,
              cart->g_rows);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int dims[2] = {cart->g_rows, cart->g_cols};
  int periods[2] = {0, 0};
  int reorder = 0;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart->cart_comm);
  if (cart->cart_comm == MPI_COMM_NULL) {
    if (cart->my_rank == 0)
      fprintf(stderr, "MPI_Cart_create failed\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Comm_rank(cart->cart_comm, &cart->cart_rank);

  int coords[2];
  MPI_Cart_coords(cart->cart_comm, cart->cart_rank, 2, coords);
  cart->p_row = coords[0];
  cart->p_col = coords[1];

  int keep[2] = {0, 1};
  MPI_Cart_sub(cart->cart_comm, keep, &cart->row_comm);

  keep[0] = 1;
  keep[1] = 0;
  MPI_Cart_sub(cart->cart_comm, keep, &cart->col_comm);

  MPI_Comm_rank(cart->row_comm, &cart->row_rank);
}

int blocked_floyd_warshall_p_apsp(Matrix w, Matrix *buf, int b) {
  CartInfo cart;
  setup_cart(&cart, w.n);

  int n = w.n;
  int p = cart.g_rows;
  int localN = n / p;

  if (b != localN) {
    if (cart.my_rank == 0) {
      fprintf(stderr,
              "For this parallel version set b == n/sqrt(P) (got b=%d, need "
              "%d)\n",
              b, localN);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  copy_matrix(w, buf);
  if ((*buf).values == NULL) {
    fprintf(stderr, "Failed to copy initial matrix\n");
    MPI_Comm_free(&cart.row_comm);
    MPI_Comm_free(&cart.col_comm);
    MPI_Comm_free(&cart.cart_comm);
    return 1;
  }

  int B = n / b;
  int **wkk = NULL, **wkj = NULL, **wik = NULL, **wij = NULL;

  for (int k = 0; k < B; k++) {
    // dependant phase
    wkk = get_block((*buf).values, k, k, b, n);
    floyd((*buf).values, wkk, wkk, wkk, b, n);
    free(wkk);
    wkk = NULL; // <-- free immediately

    // partially dependant phase
    for (int j = 0; j < B; j++) {
      if (j == k)
        continue;
      wkj = get_block((*buf).values, k, j, b, n);
      wkk = get_block((*buf).values, k, k, b, n);
      floyd((*buf).values, wkj, wkk, wkj, b, n);
      free(wkj);
      wkj = NULL;
      free(wkk);
      wkk = NULL;
    }

    for (int i = 0; i < B; i++) {
      if (i == k)
        continue;
      wik = get_block((*buf).values, i, k, b, n);
      wkk = get_block((*buf).values, k, k, b, n);
      floyd((*buf).values, wik, wik, wkk, b, n);
      free(wkk);
      wkk = NULL;

      // independant phase
      for (int j = 0; j < B; j++) {
        if (j == k)
          continue;
        wkj = get_block((*buf).values, k, j, b, n);
        wij = get_block((*buf).values, i, j, b, n);
        floyd((*buf).values, wij, wik, wkj, b, n);
        free(wkj);
        wkj = NULL;
        free(wij);
        wij = NULL;
      }
      free(wik);
      wik = NULL;
    }
  }

  MPI_Comm_free(&cart.row_comm);
  MPI_Comm_free(&cart.col_comm);
  MPI_Comm_free(&cart.cart_comm);
  return 0;
}
