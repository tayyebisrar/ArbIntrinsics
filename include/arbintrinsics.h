#ifndef ARBINTRINSICS_H
#define ARBINTRINSICS_H

#include <stdint.h>  // needed for uint64_t, uint8_t, etc.
#include <stddef.h>
#include <stdbool.h> // for bool (sign)

typedef uint64_t limb_t;

typedef enum{
    ARBINTR_TYPE_NORMAL = 0,
    ARBINTR_TYPE_ZERO,
    ARBINTR_TYPE_INF,
    ARBINTR_TYPE_NAN
} arbfloat_type_t;

typedef struct {
    limb_t *limbs;       // pointer to dynamically allocated array
    size_t length;       // how many limbs are used
    size_t capacity;     // how many limbs are allocated
    int exponent;        // unbiased exponent
    bool sign;           // 0 for +, 1 for -
    arbfloat_type_t type; // type of the arbfloat (e.g. normal, subnormal, zero, infinity, NaN)
    // int user_prec_bits; // user-defined precision bits (not a multiple of 64)
} arbfloat_t;

typedef enum{
    RNDN=0,// Round to nearest (tie to even)
    RNDZ,// Round to zero
    RNDA,// Round away from zero
    RNDU,// Round up (+inf)
    RNDD // Round down (-inf)
} roundmode;

#define ARBINTR_EXP_MAX 1024
#define ARBINTR_EXP_MIN -1024

#define ARBINTR_IS_NORMAL(x) ((x)->type == ARBINTR_TYPE_NORMAL)
#define ARBINTR_IS_INF(x) ((x)->type == ARBINTR_TYPE_INF)
#define ARBINTR_IS_NAN(x) ((x)->type == ARBINTR_TYPE_NAN)
#define ARBINTR_IS_ZERO(x) ((x)->type == ARBINTR_TYPE_ZERO)

#define ARBINTR_EXP_IS_IN_RANGE(x) \
    ((x)->exponent >= ARBINTR_EXP_MIN && (x)->exponent <= ARBINTR_EXP_MAX)

#define ARBINTR_IS_NORMALISED(x) ((x)->limbs[0] >> 63 & 1) // MSB is 1

#endif