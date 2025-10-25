#ifndef COMM_H
#define COMM_H

#include "types.h"

void setup_cart(CartInfo *cart, int n);

void gather_matrix(Matrix *local_buf, Matrix *global_buf, int b,
                   CartInfo *cart);
void distribute_matrix(Matrix w, Matrix *local_buf, int b, CartInfo *cart);

#endif // COMM_H
