#ifndef HELPERS_H
#define HELPERS_H
#include "arbintrinsics.h"

void normalise(arbfloat_t *x);
void left_shift(arbfloat_t *x, int n);
void right_shift(arbfloat_t *x, int n, uint8_t *grt);

void copy_arbfloat(arbfloat_t* dest, const arbfloat_t *src);
void append_limb(arbfloat_t *x, limb_t limb);
int compare_arbfloat_t(arbfloat_t *x, arbfloat_t *y);

void check_overflow(arbfloat_t *x, roundmode roundmode);

#endif