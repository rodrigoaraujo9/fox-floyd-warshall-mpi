#ifndef PARALLEL_FW_H
#define PARALLEL_FW_H

#include "types.h"

void setup_cart(CartInfo *cart, int n);

int send_block_p(Matrix m, int n, CartInfo cart);

int blocked_floyd_warshall_p_apsp(Matrix w, Matrix *buf, int b);

#endif // PARALLEL_FW_H
