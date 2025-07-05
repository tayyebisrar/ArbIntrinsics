#ifndef TESTHELPERS_H
#define TESTHELPERS_H
#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <math.h>

void generate_random_arbfloat(arbfloat_t *x, size_t length);
void display_arbfloat(arbfloat_t *x);
void tobinary(limb_t *n);
limb_t rand64(void);

#endif