#ifndef OPERATIONS_H
#define OPERATIONS_H
#include "arbintrinsics.h"

void addnaive(arbfloat_t *x, arbfloat_t *y,arbfloat_t *res, roundmode roundmode);
void multnaive(arbfloat_t *x, arbfloat_t *y, arbfloat_t *res, roundmode roundmode);
void multkaratsuba(arbfloat_t *x, arbfloat_t *y,arbfloat_t *res);
// void multtoomcook(arbfloat_t *x, arbfloat_t *y, int k, arbfloat_t *z);

#endif