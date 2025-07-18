#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <math.h>
#include <string.h>
#include "fphelpers.h"
#include "testhelpers.h"

#define SIMD_ALL_ZEROES _mm256_setzero_si256()

void add_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov);

void addnaive(arbfloat_t *x, arbfloat_t *y,arbfloat_t *res){
    // x and y themselves are unchanged, use copies when shifting
    arbfloat_t *op1 = (x->exponent >= y->exponent) ? x : y; // Useful for easy swapping
    arbfloat_t *op2 = (x->exponent >= y->exponent) ? y : x;
    arbfloat_t *op2_copy;

    res->exponent = op1->exponent; // set the exponent of the result
    res->sign = op1->sign;
    res->length = op1->length;
    res->capacity = op1->capacity;
    res->limbs = malloc(res->capacity * sizeof(limb_t));
    uint8_t grt;

    int d = op1->exponent - op2->exponent; // exp(op1) is always >= exp(op2), so d >= 0
    if (d > 0) { // if difference in exps, align op2 by right shifting
        // use a (deep) copy of op2 instead to avoid changing op2, so reassign the op2 pointer to a copy
        printf("Before shift and reassign op2 to copy:\n");
        printf("op1:");
        display_arbfloat_limbs(op1, 1);
        printf("op2:");
        display_arbfloat_limbs(op2, 1);
        op2_copy = copy_arbfloat(op2);
        op2 = op2_copy; // now op2 points to the copy
        right_shift(op2, d, &grt);
        printf("After shift and reassign op2 to copy:\n");
        printf("op1:");
        display_arbfloat_limbs(op1, 1);
        printf("op2:");
        display_arbfloat_limbs(op2, 1);
    }

    // if signs are same, add limbs. If signs different, subtract
    if (op1->sign ^ op2->sign) {
        //sub_limbwise(&op1, &op2, &z);
    } else {
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
    }
    if (d > 0) {
        free(op2_copy->limbs);
        free(op2_copy);
    }
    // normalise the result
    //normalise(res);
}

void add_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov){
    bool overflow = false, temp = false; // indicates if need to carry something for next add. At start there is no previous carry
    for (int i=op1->length-1; i >= 0; i--){ // for the last 4 or less limbs
        printf("Adding limbs %d\n", i);
        limb_t result = op1->limbs[i] + op2->limbs[i];
        printf("op1: %llu\n+\nop2: %llu\n=res: %llu\n", op1->limbs[i], op2->limbs[i], result);
        temp = (result < op1->limbs[i]);
        result += overflow;
        overflow = (temp || result < op1->limbs[i]);
        if (overflow) {
            printf("Overflow\n");
        }
        else{
            printf("No Overflow\n");
        }
        printf("POST OVERFLOW res:%llu\n\n", result);
        res->limbs[i] = result;
    }
    *ov = overflow;
}

void sub_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov){

}

/*

        __m256i carryinfo = _mm256_cmpgt_epi64(op1curr4, oresult); // If any overflow, the mirror limb in carryinfo will be all 1s
        simdtemp = carryinfo'
        oresult = _mm256_add_epi64(oresult, overflow); // op1 + op2 + previous set's overflow (only for 4th limb)
        __m256i carryinfo = _mm256_or_si256(carryinfo, _mm256_cmpgt_epi64(op1curr4, oresult)); // equivalent to overflow = (temp || result < op1->limbs[i]);

        // right shift carry info to get 1 if overflow, 0 otherwise. 
        // then left shift + wraparound the limbs so they're above the correct limb for adding
        carryinfo = _mm256_srli_epi64(carryinfo, 63); // overflow = 1, none = 0
        __m256i permuted = _mm256_permute4x64_epi64(carryinfo, 0x39);

        __m256i simdtemp = _mm256_extract_epi64(permuted, 0x0); // hold onto the lastlimb for next set of 4

        //carryinfo = _mm256_and_si256(carryinfo, carrymask); // clear the last lane
        oresult = _mm256_add_epi64(oresult, carryinfo); // add overflows
        
        __m256i carryinfo = _mm256_cmpgt_epi64(op1curr4, oresult); // if new sum with overflow, overflows

        //carryinfo = _mm256_and_si256(carryinfo, (111111....64x0s)); // need to blend carryinfo with 0s to remove the LSL overflow
        
        // if it WAS overflow, need to remember for the next set 

        // overflow = (_mm256_extract_epi64(carryinfo, 0) != 0); // If next set of limbs needs overflow, last (0th) limb = 1
        
        */