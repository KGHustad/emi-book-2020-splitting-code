#ifndef _MATRIX_CSR_H
#define _MATRIX_CSR_H

#include <stdint.h>
#include <stdio.h>

typedef int32_t idx_t;

struct matrix_csr
{
    uint64_t num_rows;
    uint64_t num_cols;
    uint64_t num_nonzeros;
    idx_t *row_ptr;
    idx_t *col;
    double *val;
};


void read_csr_binary(struct matrix_csr *m, FILE* f);
void sparse_mvm_csr(struct matrix_csr *m, double *B, double *C);

#endif
