#include <stdio.h>
#include <stdlib.h>
#include "arbintrinsics.h"
#include "operations.h"
#include "fphelpers.h"

#define ULP(x) x->limbs[x->length - 1] & 1

void add_rnd(arbfloat_t *x, bool rnd){
    bool overflow = rnd;
    int end= ARBINTR_LENGTH(x);
    int i = 0;
    do{
        limb_t res = x->limbs[i] += overflow;
        overflow = 0;
        if (res < x->limbs[i]){ // overflow again
            overflow = 1;
        }
        x->limbs[i] = res;
    }while (overflow && i<end);
    if (overflow){
        right_shift(x, 1, NULL);
        x->limbs[end-1] += (1ULL << 63);
    }
    return; // ditch any new GRT bits
}

void roundup(arbfloat_t *x, uint8_t grt){
    /* ROUND TO POSITIVE INF
    If any of GRT is set:
    - If x is positive, round up
    - If x is negative, don't change */

    if (grt){ // grt != 000
        if (x->sign){ // negative
            return; // don't change
        }
        else{
            add_rnd(x, 1);
        }
    }
    return;
}

void rounddown(arbfloat_t *x, uint8_t grt){
    /* ROUND TO NEGATIVE INF
    If any of GRT is set:
    - If x is positive, don't change 
    - If x is negative, round 'up' (add 1 to mantissa)*/
    
    if (grt) {
        if (x->sign){
            add_rnd(x, 1);
        }
        else{
            return;
        }
    }
    return;
}

void roundtonearest(arbfloat_t *x, uint8_t grt){
    /*ROUND TO NEAREST - TIE TO EVEN
    Rounds up if x is half a ulp extra or more,
    i.e. rnd = G(R + T) for rounding up, and
    rnd = G(R+T)'L for rounding to even, since if L=1,
    needs to become L=0. Put together,
    rnd = G(R + T + L)*/

    bool rnd = (((grt << 1) | (uint8_t)ULP(x)) > 8); // GRT becomes GRTL, so GRT > 0b1000 (8) <=> G(R+T+L)

    if (rnd){
        add_rnd(x, 1);
    }
}

void roundtozero(arbfloat_t *x, uint8_t grt){
    /*ROUND TO ZERO - TRUNCATE*/
    return;
}

void roundawayfromzero(arbfloat_t *x, uint8_t grt){
    /*ROUND AWAY FROM ZERO
    If x is positive, round to infinity (up)
    If x is negative, round to neg. infinity (down) */;
    if (x->sign){
        rounddown(x, grt);
    }
    else{
        roundup(x, grt);
    }
}

void round_arbfloat(arbfloat_t *x, uint8_t grt, roundmode mode){
    // Rounds the arbfloat x according to the specified rounding mode
    if (!ARBINTR_IS_NORMAL(x)) return;
    switch (mode) {
        case RNDN:
            roundtonearest(x, grt);
            break;
        case RNDD:
            rounddown(x, grt);
            break;
        case RNDU:
            roundup(x, grt);
            break;
        case RNDZ:
            roundtozero(x, grt);
            break;
        case RNDA:
            roundawayfromzero(x, grt);
            break;
        default:
            fprintf(stderr, "Unknown rounding mode\n");
            exit(EXIT_FAILURE);
    }
}