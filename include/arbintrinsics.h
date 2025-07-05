#ifndef ARBINTRINSICS_H
#define ARBINTRINSICS_H

#include <stdint.h>  // needed for uint64_t, uint8_t, etc.
#include <stddef.h>
#include <stdbool.h> // for bool (sign)

typedef uint64_t limb_t;

typedef struct {
    limb_t *limbs;       // pointer to dynamically allocated array
    size_t length;       // how many limbs are used
    size_t capacity;     // how many limbs are allocated
    int exponent;        // unbiased exponent
    bool sign;            // 0 for +, 1 for -
} arbfloat_t;

typedef enum{
    RNDN=0,// Round to nearest (tie to even)
    RNDZ,// Round to zero
    RNDU,// Round up (+inf)
    RNDD // Round down (-inf)
}roundmode;

#endif