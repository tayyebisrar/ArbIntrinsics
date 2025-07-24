#ifndef HELPERS_H
#define HELPERS_H
#include "arbintrinsics.h"

void roundtonearest(arbfloat_t *x, uint8_t grt);
void roundtozero(arbfloat_t *x, uint8_t grt);
void roundawayfromzero(arbfloat_t *x, uint8_t grt);
void roundup(arbfloat_t *x, uint8_t grt);
void rounddown(arbfloat_t *x, uint8_t grt);

#endif