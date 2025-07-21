#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <math.h>
#include <string.h>
#include "fphelpers.h"

#define SIMD_ALL_ZEROES _mm256_setzero_si256()

void add_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov);
void sub_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res);

void addnaive(arbfloat_t *x, arbfloat_t *y,arbfloat_t *res, uint8_t* grtp){
    // TODO: Rounding
    // TODO: Change normalise to happen after rounding too (?)
    // x and y themselves are unchanged, use copies when shifting
    arbfloat_t *op1 = (x->exponent >= y->exponent) ? x : y; // Useful for easy swapping
    arbfloat_t *op2 = (x->exponent >= y->exponent) ? y : x;
    arbfloat_t *op2_copy, *op1_copy;

    res->exponent = op1->exponent; // set the exponent of the result
    res->sign = op1->sign;
    res->length = op1->length;
    res->capacity = op1->capacity;
    res->limbs = malloc(res->capacity * sizeof(limb_t));
    uint8_t grt_out = 0;

    int d = op1->exponent - op2->exponent; // exp(op1) is always >= exp(op2), so d >= 0
    
    // if signs are same, add limbs. If signs are different, subtract
    // we always take the sign of the larger magnitude if different

    if (x->sign ^ y->sign){
        uint8_t grt = 0;
        if (op1->exponent != op2->exponent){ // since exponents are different, and op1's is larger, it can always be assumed here op1 > op2
            op2_copy = copy_arbfloat(op2);
            op2 = op2_copy; // now op2 points to the copy
            op1_copy = copy_arbfloat(op1);
            op1 = op1_copy;
            right_shift(op2, d, &grt); // using an aligned subtrahend
            append_limb(op1, 0); // needed to accurately perform subtraction with GRT
            append_limb(op2, ((limb_t)(grt&0x07)<<61)); // places GRT at the left side of the new limb, i.e. GRT0000....
            append_limb(res, 0);
            sub_limbwise(op1, op2, res); // larger - smaller
            
            normalise(res); // normalise result
            // recompute GRT bits from extra limb
            uint8_t tgrt = 0;
            limb_t extralimb = res->limbs[res->length-1];
            tgrt |= ((extralimb >> 63) & 1) << 2; // new G bit
            tgrt |= ((extralimb >> 62) & 1) << 1; // new R bit
            tgrt |= (extralimb & ((1ULL << (62)) - 1)) ? 1 : 0;

            // remove the extra limbs
            op1->length--;
            op2->length--;
            res->length--;

            grt_out = tgrt;
        }
        else{ // they are of equal exponent - all is good, and no extra limbs are needed since they're the same exponent (no grt)
            if (compare_arbfloat_t(op1, op2) == -1){  // if op1 is actually < op2 (based on mantissas)
                sub_limbwise(op2, op1, res);
                res->sign = op2->sign; // the sign should be changed to the larger in this case (op2)
            }
            else{
                sub_limbwise(op1, op2, res); 
            }
            grt_out = grt; // no grt, so will be 0
        }
    }
    else {
        uint8_t grt = 0;
        if (d > 0) { // if difference in exps, align op2 by right shifting
            // use a (deep) copy of op2 instead to avoid changing op2, so reassign the op2 pointer to a copy
            op2_copy = copy_arbfloat(op2);
            op2 = op2_copy; // now op2 points to the copy
            right_shift(op2, d, &grt);
        }
        bool overflowed;
        add_limbwise(op1, op2, res, &overflowed);
        uint8_t tgrt = 0;
        if (overflowed){
            right_shift(res, 1, &tgrt); // first limb will go bbbbbbbbbb to 0bbbbbbbbbb
            // TGRT bits = old grt bits' last 2 values bitwise ORed and placed in 2nd-to-last position
            // then right shifted once, then bitwise OR'd with tgrt
            // i.e GRT -> 0G(R+T) -> TG(R+T) -> new GRT. 
            grt = (grt >> 1) | (grt & 1) | tgrt;
            res->limbs[0] += (1ULL << 63); // since we needed overflow, first limb goes 0bbbbbbbbbb to 1bbbbbbbbbb
        }
        grt_out = grt;
    }
    if (op2_copy) {
        free(op2_copy->limbs);
        free(op2_copy);
    }
    if (op1_copy) {
        free(op1_copy->limbs);
        free(op1_copy);
    }
    if (grtp) *grtp = grt_out;

    // normalise(res) here might change the GRT? 
    return;
}  

void add_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov){
    bool overflow = false, temp = false; // indicates if need to carry something for next add. At start there is no previous carry
    for (int i=op1->length-1; i >= 0; i--){ // for the last 4 or less limbs
        limb_t result = op1->limbs[i] + op2->limbs[i];
        temp = (result < op1->limbs[i]);
        result += overflow;
        overflow = (temp || result < op1->limbs[i]);
        res->limbs[i] = result;
    }
    *ov = overflow;
}

void sub_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res){
    bool borrow = false, temp = false;
    for (int i=op1->length-1; i>=0; i--){
        limb_t result = op1->limbs[i] - op2->limbs[i];
        temp = (result > op1->limbs[i]);
        result -= borrow;
        borrow = (temp || result > op1->limbs[i]);
        res->limbs[i] = result;
    }
}