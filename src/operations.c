#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <math.h>
#include "fphelpers.h"

void addnaive(arbfloat_t *x, arbfloat_t *y,arbfloat_t *z){
    // REMEMBER TO CONSIDER LATER: x and y themselves are changed, might
    // want to use copies instead
    arbfloat_t *op1 = (x->exponent >= y->exponent) ? x : y; // Useful for easy swapping
    arbfloat_t *op2 = (x->exponent >= y->exponent) ? y : x;
    
    z->exponent = op1->exponent; // set the exponent of the result
    z->sign = op1->sign;
    z->length = op1->length;
    z->capacity = op1->capacity;

    int d = op1->exponent - op2->exponent; // exp(op1) is always >= exp(op2), so d >= 0
    if (d > 0) { // if difference in exps
        // align op2 by right shifting
        uint8_t grt;
        right_shift(op2, d, &grt);
        //op2->exponent -= d; // return exponent to original - the original operand shouldn't change
    }

    // if signs are same, add limbs. If signs different, subtract
    if (z->sign == op1->sign) {
        //add_limbwise(op1, op2, z);
    } else {
        //sub_limbwise(op1, op2, z);
    }
    
    // normalise the result
    //normalise(z);
}