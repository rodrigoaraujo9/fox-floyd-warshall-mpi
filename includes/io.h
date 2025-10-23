#ifndef IO_H
#define IO_H

#include "types.h"

Matrix read_input_matrix_from_file(char *file_path);
Matrix read_output_matrix_from_file(char *file_path, int n);
void print_matrix(Matrix matrix, char *title);

#endif // IO_H
