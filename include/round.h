#ifndef ROUND_H
#define ROUND_H
#include "arbintrinsics.h"

void round_arbfloat(arbfloat_t *x, uint8_t grt, roundmode mode);
void roundtonearest(arbfloat_t *x, uint8_t grt);
void roundtozero(arbfloat_t *x, uint8_t grt);
void roundawayfromzero(arbfloat_t *x, uint8_t grt);
void roundup(arbfloat_t *x, uint8_t grt);
void rounddown(arbfloat_t *x, uint8_t grt);

#endif