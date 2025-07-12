#ifndef HELPERS_H
#define HELPERS_H
#include "arbintrinsics.h"

void normalise(arbfloat_t *x);
void left_shift(arbfloat_t *x, int n);
void right_shift(arbfloat_t *x, int n, uint8_t *grt);

#endif