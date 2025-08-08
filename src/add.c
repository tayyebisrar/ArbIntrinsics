#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <operations.h>
#include <fphelpers.h>

void arbintr_add(arbfloat_t *a, arbfloat_t *b, arbfloat_t* res, roundmode roundmode){
    // Addition itself is handled by addnaive
    // This function checks special cases, calling addnaive if both normal

    // Handle special cases - start with non-normal cases
    if (ARBINTR_IS_NAN(a) || ARBINTR_IS_NAN(b)) {
        // NaN + anything = NaN
        res->type = ARBINTR_TYPE_NAN;
        return;
    }
    else if (ARBINTR_IS_INF(a)){
        if (!ARBINTR_IS_INF(b) || a->sign == b->sign){ // Inf + +Inf, Inf + b, -Inf - Inf
            // Inf + anything = Inf (if same sign)
            res->type = ARBINTR_TYPE_INF;
            res->sign = a->sign;
            return;
        }
        else{
            // Inf - Inf
            res->type = ARBINTR_TYPE_NAN; // Inf - Inf = NaN
            return;
        }
    }
    else if (ARBINTR_IS_INF(b)){ // a + Inf, a - Inf (every Inf case is now handled)
        res->type = ARBINTR_TYPE_INF;
        res->sign = b->sign;
        return;
    }
    else if (ARBINTR_IS_ZERO(a)){
        if (ARBINTR_IS_ZERO(b)) {
            // 0 + 0 = 0, -0 only when both are zero, but if RNDN, -0 if either is 0
            res->type = ARBINTR_TYPE_ZERO;
            res->sign = (roundmode == RNDN) ? (a->sign || b->sign) : (a->sign && b->sign);
            return;
        }
        else {
            // 0 + b = b
            copy_arbfloat(res, b, NULL); // res = b
            return;
        }
    }
    else if (ARBINTR_IS_ZERO(b)){
        // a + 0 = a
        copy_arbfloat(res, a, NULL);
        return;
    }
    else{
        // Make sure both are normalised
        if (!ARBINTR_IS_NORMALISED(a) || !ARBINTR_IS_NORMALISED(b)) {
            normalise(a);
            check_overflow(a);
            normalise(b);
            check_overflow(b);
            arbintr_add(a, b, res, roundmode); // call again with normalised operands
            return;
        }
        addnaive(a, b, res, roundmode); // both normal, call addnaive
        res->type = ARBINTR_TYPE_NORMAL; // set type to normal
        check_overflow(res); // update to special value if overflow/underflow
        return;
    }
}