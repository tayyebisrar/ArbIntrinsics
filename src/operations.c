#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <math.h>
#include "fphelpers.h"

#define SIMD_ALL_ZEROES _mm256_setzero_si256()

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
        //add_limbwise(&op1, &op2, &z);
    } else {
        //sub_limbwise(&op1, &op2, &z);
    }
    
    // normalise the result
    //normalise(z);
}

void add_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* z, bool *ov){
    __m256i zeroes = SIMD_ALL_ZEROES;
    for (int i=op1->length-4; i>3; i--){
        __m256i op1curr4 = _mm256_loadu_si256((const __m256i *)(op1->limbs + i));
        __m256i op2curr4 = _mm256_loadu_si256((const __m256i *)(op2->limbs + i));

        __m256i result = _mm256_add_epi64(op1curr4, op2curr4);

        __m256i carryinfo = _mm256_cmpgt_epi64(op1curr4, result);
        carryinfo = _mm256_or_si256(carryinfo, _mm256_cmpgt_epi64(op2curr4, result)); // if (result < op1->limbs[i] || result < op2->limbs[i])

        carryinfo = _mm256_srli_epi64(carryinfo, 63); // overflow = 1, none = 0

        carryinfo = _mm256_permute4x64_epi64(carryinfo, 0x39); // left shift overflows over the correct limbs for adding
        carryinfo = _mm256_and_si256(carryinfo, ()); // need to blend carryinfo with 0s to remove the LSL overflow
        
        // if it WAS overflow, need to remember for the next set 


        
    }

    bool overflow = false;
    for (int i=(op1->length < 4) ? op1->length-1 : 3 + op1->length % 4; i >= 0; i--){
        limb_t result = op1->limbs[i] + op2->limbs[i] + overflow;
        overflow = false;
        if (result < op1->limbs[i] || result < op2->limbs[i]){
            overflow = true;
        }
        z->limbs[i] = result;
    }
    *ov = overflow;
}

/*
bool temp = false
for (int i=(op1->length >= 4) ? 4 :op1->length; i >= 0; i--){
        limb_t result = op1->limbs[i] + op2->limbs[i];
        if (result < op1->limbs[i]){
            temp = true;
        }
        result += overflow; // cant overflow twice, so this is fine
        overflow = temp;
        z->limbs[i] = result;
    }

    *ov = overflow;
*/