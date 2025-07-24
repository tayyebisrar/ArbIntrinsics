#include <stdio.h>
#include <stdlib.h>
#include "arbintrinsics.h"
#include "operations.h"
#include "fphelpers.h"

#define ULP(x) x->limbs[x->length - 1] & 1

void add_rnd(arbfloat_t *x, bool rnd){
    bool overflow = rnd;
    int i=x->length-1;
    do{
        limb_t res = x->limbs[i] += overflow;
        overflow = 0;
        if (res < x->limbs[i]){ // overflow again
            overflow = 1;
        }
        x->limbs[i] = res;
    }while (overflow && i>=0);
    if (overflow){
        right_shift(x, 1, NULL);
        x->limbs[0] += (1ULL << 63);
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

    bool rnd = (((grt << 1) && ULP(x)) > 8); // GRT becomes GRTL, so GRT > 0b1000 (8) <=> G(R+T+L)

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

