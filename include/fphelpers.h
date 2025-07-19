#ifndef HELPERS_H
#define HELPERS_H
#include "arbintrinsics.h"

void normalise(arbfloat_t *x);
void left_shift(arbfloat_t *x, int n);
void right_shift(arbfloat_t *x, int n, uint8_t *grt);

arbfloat_t* copy_arbfloat(const arbfloat_t *src);
void append_limb(arbfloat_t *x, limb_t limb);
#endif