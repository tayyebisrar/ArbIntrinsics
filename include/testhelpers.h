#ifndef TESTHELPERS_H
#define TESTHELPERS_H
#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <math.h>

void generate_random_arbfloat(arbfloat_t *x, size_t length);
void display_arbfloat_limbs(arbfloat_t *x, int oneline);
void display_normalised(arbfloat_t *x);
void tobinary(limb_t n);
limb_t rand64(void);
void print_m256i_binary(__m256i val, int gaps);
void set_arbfloat(arbfloat_t *x, limb_t *limbs, size_t l, int exponent, bool sign);

#endif