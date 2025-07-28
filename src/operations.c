#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <math.h>
#include <string.h>
#include "fphelpers.h"
#include "testhelpers.h"
#include "round.h"

#define SIMD_ALL_ZEROES _mm256_setzero_si256()

void add_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov);
void sub_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res);
void add_simd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov);
void sub_simd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res);

void addnaive(arbfloat_t *x, arbfloat_t *y,arbfloat_t *res, roundmode roundmode){
    // Method: Align -> compute grt -> add/subtract -> normalise -> round
    // Assumes x and y are normalised, and are both normal.
    // TODO: Find an unsigned version of _mm256_cmpgt_epi64 to correctly compare limbs. Could use XORs or AVX512
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
            op2_copy = malloc(sizeof(arbfloat_t));
            copy_arbfloat(op2_copy, op2);
            op2 = op2_copy; // now op2 points to the copy
            op1_copy = malloc(sizeof(arbfloat_t));
            copy_arbfloat(op1_copy, op1);
            op1 = op1_copy;
            right_shift(op2, d, &grt); // using an aligned subtrahend
            uint8_t subtracted_grt = (-grt + 1) & 0x07; // 00000111 = 0x07, computes 000 - GRT
            sub_limbwise(op1, op2, res); // larger - smaller
            append_limb(res, ((limb_t)(subtracted_grt&0x07)<<61)); // places 000-GRT in the last limb of result
            normalise(res); // normalise result

            // recompute GRT bits from extra limb
            uint8_t tgrt = 0;
            limb_t extralimb = res->limbs[res->length-1];
            tgrt |= ((extralimb >> 63) & 1) << 2; // new G bit
            tgrt |= ((extralimb >> 62) & 1) << 1; // new R bit
            tgrt |= (extralimb & ((1ULL << (62)) - 1)) ? 1 : 0;

            // remove the extra limbs
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
            op2_copy = malloc(sizeof(arbfloat_t));
            copy_arbfloat(op2_copy, op2);
            op2 = op2_copy; // now op2 points to the copy
            right_shift(op2, d, &grt);
        }
        bool overflowed;
        add_limbwise(op1, op2, res, &overflowed);
        if (overflowed){
            uint8_t tgrt = 0; // will either be 0 or 1, based on ulp
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
        if (op2_copy->limbs) {
            free(op2_copy->limbs); // free the limbs of the copy
        }
        free(op2_copy);
    }
    if (op1_copy) {
        if (op1_copy->limbs) {
            free(op1_copy->limbs);
        }
        free(op1_copy);
    }

    if (roundmode){
        printf("Before rounding: ");
        display_arbfloat_limbs(res, 1);
        round_arbfloat(&res, grt_out, roundmode);
        printf("After rounding: ");
        display_arbfloat_limbs(res, 1);
    }
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

void add_simd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov){
    __m256i mask = _mm256_set_epi64x(0, -1, -1, -1);
    __m256i next_carry_mask = _mm256_set_epi64x(-1, 0, 0, 0);
    __m256i new_result;
    __m256i carryover = SIMD_ALL_ZEROES; // starts with no carry
    for (int i = op1->length - 4; i >= 0; i -= 4) {
        __m256i a = _mm256_loadu_si256((const __m256i *)(op1->limbs + i));
        __m256i b = _mm256_loadu_si256((const __m256i *)(op2->limbs + i));
        __m256i result = _mm256_add_epi64(a, b); // add the 4 limbs
        __m256i cmp_gt = _mm256_cmpgt_epi64(a, result); // <=> result < a
        result = _mm256_add_epi64(result, carryover); // add the carry
        __m256i cmp_gt_with_carry = _mm256_or_si256(cmp_gt, _mm256_cmpgt_epi64(a, result));
        cmp_gt_with_carry = _mm256_permute4x64_epi64(cmp_gt_with_carry, 0x39); // moves each limb left by 1, so all limbs (except the first) can be directly added
        cmp_gt_with_carry = _mm256_srli_epi64(cmp_gt_with_carry, 63); // shifts 0xFFF.. to 0x000..1, so any carry limbs can now be added
        __m256i maskedcarry = _mm256_and_si256(cmp_gt_with_carry, mask);
        carryover = _mm256_and_si256(cmp_gt_with_carry, next_carry_mask); // carryover is the carry for the next iteration);
        while (!_mm256_testz_si256(maskedcarry, maskedcarry)){ // while there is a carry, repeat. this instruction returns 1 if carry is zero, 0 otherwise
            // Propagate the carries from within the same set of 4 limbs
            new_result = _mm256_add_epi64(result, maskedcarry); // each limb except limb 0 + either 0 or 1
            cmp_gt = _mm256_cmpgt_epi64(result, new_result); // if *another* overflow happens, newresult.limb < result.limb
            result = new_result;
            cmp_gt = _mm256_permute4x64_epi64(cmp_gt, 0x39); // left shift
            cmp_gt = _mm256_srli_epi64(cmp_gt, 63); // shift to the right to get 0xFFF.. to 0x000..1
            maskedcarry = _mm256_and_si256(cmp_gt, mask);
            carryover = _mm256_or_si256(carryover, _mm256_and_si256(cmp_gt, next_carry_mask)); // limb0 out of 0-3 for the next iteration
        }
        _mm256_storeu_si256((__m256i *)(res->limbs + i), result); // store the result
    }

    // If carryover is not zero, there is an overall overflow
    bool overflow = false, temp = false; 
    overflow = !_mm256_testz_si256(carryover, carryover);
    // Handle any remaining limbs (less than 4)
    for (int i=(op1->length % 4) -1; i >= 0; i--){
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

void sub_simd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res){
    __m256i mask = _mm256_set_epi64x(0, -1, -1, -1);
    __m256i next_borrow_mask = _mm256_set_epi64x(-1, 0, 0, 0);
    __m256i new_result;
    __m256i borrowover = SIMD_ALL_ZEROES; // starts with no borrow
    for (int i = op1->length - 4; i >= 0; i -= 4) {
        __m256i a = _mm256_loadu_si256((const __m256i *)(op1->limbs + i));
        __m256i b = _mm256_loadu_si256((const __m256i *)(op2->limbs + i));
        __m256i result = _mm256_sub_epi64(a, b); // subtract the 4 limbs
        __m256i cmp_gt = _mm256_cmpgt_epi64(result, a); // <=> result > a
        result = _mm256_sub_epi64(result, borrowover); // subtract the borrow
        __m256i cmp_gt_with_borrow = _mm256_or_si256(cmp_gt, _mm256_cmpgt_epi64(result, a));
        cmp_gt_with_borrow = _mm256_permute4x64_epi64(cmp_gt_with_borrow, 0x39); // moves each limb left by 1, so all limbs (except the first) can be directly added
        cmp_gt_with_borrow = _mm256_srli_epi64(cmp_gt_with_borrow, 63); // shifts 0xFFF.. to 0x000..1, so any borrow limbs can now be added
        __m256i maskedborrow = _mm256_and_si256(cmp_gt_with_borrow, mask);
        borrowover = _mm256_and_si256(cmp_gt_with_borrow, next_borrow_mask); // borrowover is the borrow for the next iteration
        while (!_mm256_testz_si256(maskedborrow, maskedborrow)){ // while there is a borrow, repeat. this instruction returns 1 if borrow is zero, 0 otherwise
        // Propagate the borrows from within the same set of 4 limbs
            new_result = _mm256_sub_epi64(result, maskedborrow); // each limb except limb 0 + either 0 or 1
            cmp_gt = _mm256_cmpgt_epi64(new_result, result); // if *another* overflow happens, newresult.limb > result.limb
            result = new_result;
            cmp_gt = _mm256_permute4x64_epi64(cmp_gt, 0x39); // left shift
            cmp_gt = _mm256_srli_epi64(cmp_gt, 63); // shift to the right to get 0xFFF.. to 0x000..1
            maskedborrow = _mm256_and_si256(cmp_gt, mask);
            borrowover = _mm256_or_si256(borrowover, _mm256_and_si256(cmp_gt, next_borrow_mask)); // limb0 out of 0-3 for the next iteration
        }
        _mm256_storeu_si256((__m256i *)(res->limbs + i), result); // store the result
    }

    // If borrowover is not zero, there is an overall overflow
    bool borrow = false, temp = false;
    borrow = !_mm256_testz_si256(borrowover, borrowover);
    // Handle any remaining limbs (less than 4)
    for (int i=(op1->length % 4) -1; i >= 0; i--){
        limb_t result = op1->limbs[i] - op2->limbs[i];
        temp = (result > op1->limbs[i]);
        result -= borrow;
        borrow = (temp || result > op1->limbs[i]);
        res->limbs[i] = result;
    }
}