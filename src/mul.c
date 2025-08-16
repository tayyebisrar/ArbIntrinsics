#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <operations.h>
#include <fphelpers.h>

void arbintr_mul(arbfloat_t *a, arbfloat_t *b, arbfloat_t* res, roundmode roundmode){
    // Addition itself is handled by addnaive
    // This function checks special cases, calling addnaive if both normal

    // Handle special cases - start with non-normal cases
    if (ARBINTR_IS_NAN(a) || ARBINTR_IS_NAN(b)) {
        // NaN * anything = NaN
        res->type = ARBINTR_TYPE_NAN;
        return;
    }
    else if (ARBINTR_IS_INF(a) || ARBINTR_IS_INF(b)){
        // Inf * anything = Inf
        res->type = ARBINTR_TYPE_INF;
        res->sign = a->sign ^ b->sign;
    }
    else if (ARBINTR_IS_ZERO(a) || ARBINTR_IS_ZERO(b)){
        // 0 * anything = Zero
        res->type = ARBINTR_TYPE_ZERO;
    }
    else{
        // Make sure both are normalised
        if (!ARBINTR_IS_NORMALISED(a) || !ARBINTR_IS_NORMALISED(b)) {
            normalise(a);
            check_overflow(a);
            normalise(b);
            check_overflow(b);
            arbintr_mul(a, b, res, roundmode); // call again with normalised operands
            return;
        }
        multnaive(a, b, res, roundmode); // both normal, call addnaive
        res->type = ARBINTR_TYPE_NORMAL; // set type to normal
        check_overflow(res); // update to special value if overflow/underflow
        return;
    }
}