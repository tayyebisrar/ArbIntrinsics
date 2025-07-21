#ifndef OPERATIONS_H
#define OPERATIONS_H
#include "arbintrinsics.h"

void addnaive(arbfloat_t *x, arbfloat_t *y,arbfloat_t *res, uint8_t* grtp);
void multnaive(arbfloat_t *x, arbfloat_t *y, arbfloat_t *z);
void multkaratsuba(arbfloat_t *x, arbfloat_t *y,arbfloat_t *z);
void multtoomcook(arbfloat_t *x, arbfloat_t *y, int k, arbfloat_t *z);

#endif