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
#define KARATSUBA_THRESHOLD 8

void add_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov);
void sub_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool borrowin);
void multlimbwise(arbfloat_t *op1, arbfloat_t* op2, arbfloat_t* res);
void add_simd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov);
void sub_simd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool borrowin);
void multsimd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t *res);

static __m256i SIGNBIT_MASK = (__m256i){
    0x8000000000000000ULL,
    0x8000000000000000ULL,
    0x8000000000000000ULL,
    0x8000000000000000ULL
};

static inline __m256i _mm256_cmpgt_epu64(__m256i a, __m256i b) {
    // Unsigned version of _mm256_cmpgt_epi64
    return _mm256_cmpgt_epi64(_mm256_xor_si256(a, SIGNBIT_MASK), _mm256_xor_si256(b, SIGNBIT_MASK));
}

void addnaive(arbfloat_t *x, arbfloat_t *y,arbfloat_t *res, roundmode roundmode){
    // Method: Align -> compute grt -> add/subtract -> normalise -> round
    // Assumes x and y are normalised, and are both normal.
    //__m512i g = _mm512_set_epi64(0, 0, 0, 0, 0, 0, 0, 0);
    arbfloat_t *op1 = (x->exponent >= y->exponent) ? x : y; // Useful for easy swapping
    arbfloat_t *op2 = (x->exponent >= y->exponent) ? y : x;
    arbfloat_t *op2_copy = NULL, *op1_copy = NULL;

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
            uint8_t subtracted_grt = (8 - grt); // 00000111 = 0x07, computes (1)000 - GRT
            sub_limbwise(op1, op2, res, ((grt)?true:false)); // larger - smaller
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
            int gt = compare_arbfloat_t(op1, op2);
            if (gt == -1){  // if op1 is actually < op2 (based on mantissas)
                sub_limbwise(op2, op1, res, false); // false because no borrow was taken for grt
                res->sign = op2->sign; // the sign should be changed to the larger in this case (op2)
                normalise(res); // can safely normalise
            }
            else if (gt == 1){
                sub_limbwise(op1, op2, res, false);
                normalise(res); // can safely normalise
            }
            else{ // op1 is exactly equal to op2, gt == 0, res == 0
                res->type = ARBINTR_TYPE_ZERO;
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

    if (grt_out) round_arbfloat(res, grt_out, roundmode); // Round if necessary

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
        __m256i cmp_gt = _mm256_cmpgt_epu64(a, result); // <=> result < a
        result = _mm256_add_epi64(result, carryover); // add the carry
        __m256i cmp_gt_with_carry = _mm256_or_si256(cmp_gt, _mm256_cmpgt_epu64(a, result));
        cmp_gt_with_carry = _mm256_permute4x64_epi64(cmp_gt_with_carry, 0x39); // moves each limb left by 1, so all limbs (except the first) can be directly added
        cmp_gt_with_carry = _mm256_srli_epi64(cmp_gt_with_carry, 63); // shifts 0xFFF.. to 0x000..1, so any carry limbs can now be added
        __m256i maskedcarry = _mm256_and_si256(cmp_gt_with_carry, mask);
        carryover = _mm256_and_si256(cmp_gt_with_carry, next_carry_mask); // carryover is the carry for the next iteration);
        while (!_mm256_testz_si256(maskedcarry, maskedcarry)){ // while there is a carry, repeat. this instruction returns 1 if carry is zero, 0 otherwise
            // Propagate the carries from within the same set of 4 limbs
            new_result = _mm256_add_epi64(result, maskedcarry); // each limb except limb 0 + either 0 or 1
            cmp_gt = _mm256_cmpgt_epu64(result, new_result); // if *another* overflow happens, newresult.limb < result.limb
            result = new_result;
            cmp_gt = _mm256_permute4x64_epi64(cmp_gt, 0x39); // left shift
            cmp_gt = _mm256_srli_epi64(cmp_gt, 63); // shift to the right to get 0xFFF.. to 0x000..1
            maskedcarry = _mm256_and_si256(cmp_gt, mask);
            carryover = _mm256_or_si256(carryover, _mm256_and_si256(cmp_gt, next_carry_mask)); // limb0 out of 0-3 for the next iteration
        } // ewfoiuhfewuhfeiuhfewiuhewfiuhewfiuhfewiuhewfiuhewfiuhfewiuhfewiuhfeiuewfhiuewfhiuewfhiuewfhiuewfhewifuhewfiuhewfiuhewfiuh
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

void sub_limbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool borrowin){
    bool borrow = borrowin, temp = false;
    for (int i=op1->length-1; i>=0; i--){
        limb_t result = op1->limbs[i] - op2->limbs[i];
        temp = (result > op1->limbs[i]);
        result -= borrow;
        borrow = (temp || result > op1->limbs[i]);
        res->limbs[i] = result;
    }
}

void sub_simd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool borrowin){
    __m256i mask = _mm256_set_epi64x(0, -1, -1, -1);
    __m256i next_borrow_mask = _mm256_set_epi64x(-1, 0, 0, 0);
    __m256i new_result;
    __m256i borrowover = (borrowin) ? _mm256_set_epi64x(1, 0, 0, 0) : SIMD_ALL_ZEROES;
    for (int i = op1->length - 4; i >= 0; i -= 4) {
        __m256i a = _mm256_loadu_si256((const __m256i *)(op1->limbs + i));
        __m256i b = _mm256_loadu_si256((const __m256i *)(op2->limbs + i));
        __m256i result = _mm256_sub_epi64(a, b); // subtract the 4 limbs
        __m256i cmp_gt = _mm256_cmpgt_epu64(result, a); // <=> result > a
        result = _mm256_sub_epi64(result, borrowover); // subtract the borrow
        __m256i cmp_gt_with_borrow = _mm256_or_si256(cmp_gt, _mm256_cmpgt_epu64(result, a));
        cmp_gt_with_borrow = _mm256_permute4x64_epi64(cmp_gt_with_borrow, 0x39); // moves each limb left by 1, so all limbs (except the first) can be directly added
        cmp_gt_with_borrow = _mm256_srli_epi64(cmp_gt_with_borrow, 63); // shifts 0xFFF.. to 0x000..1, so any borrow limbs can now be added
        __m256i maskedborrow = _mm256_and_si256(cmp_gt_with_borrow, mask);
        borrowover = _mm256_and_si256(cmp_gt_with_borrow, next_borrow_mask); // borrowover is the borrow for the next iteration
        while (!_mm256_testz_si256(maskedborrow, maskedborrow)){ // while there is a borrow, repeat. this instruction returns 1 if borrow is zero, 0 otherwise
        // Propagate the borrows from within the same set of 4 limbs
            new_result = _mm256_sub_epi64(result, maskedborrow); // each limb except limb 0 + either 0 or 1
            cmp_gt = _mm256_cmpgt_epu64(new_result, result); // if *another* overflow happens, newresult.limb > result.limb
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

void multnaive(arbfloat_t *x, arbfloat_t *y, arbfloat_t *res, roundmode roundmode){
    // Add exponents -> Multiply limbs -> normalise -> round
    // x and y are the same length, just like with addition/subtraction
    // x and y are guaranteed to be normalised, result will be [1, 4)
    res->exponent = x->exponent + y->exponent; // set the exponent of the result
    res->sign = x->sign ^ y->sign;
    res->length = 2*(x->length); // Start off with 2p bits of precision
    res->capacity = 2*(x->length);
    res->limbs = calloc(res->capacity, sizeof(limb_t));
    uint8_t grt_out = 0;

    bool overflowed;
    multlimbwise(x, y, res);
    uint8_t gt = 0;
    lowerlimbsgt(res, &gt);
    res->length /= 2; // now ignore the lower p bits (might move into above function)
    if (overflowed){
        uint8_t tgt = 0; // temp g, same as with add, either yields 000 or 100
        right_shift(res, 1, &tgt); // right_shift adjusts exponent accordingly
        gt = (gt >> 1) | (gt & 1) | (tgt >> 1); // update GT to newG(oldG+oldT)
        res->limbs[0] += (1ULL << 63); // since we had overflow, first limb, 0bbbbbbbbbb to 1bbbbbbbbbb
    } // this is the 'normalise'

    round_arbfloat(res, gt, roundmode);
}


void multlimbwise(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t *res){
    // Method: Naive, schoolbook multiply
    limb_t carry = 0, sum = 0;
    for (int i=op1->length-1; i>=0; i--){
        for (int j=op2->length-1; j >= 0; j--){
            int idx = i + j + 1;
            // get the high and low parts of the multiplication
            __uint128_t res128 = ((__uint128_t)op1->limbs[i] * (__uint128_t)op2->limbs[j]); // 81
            sum = (limb_t)res128 + res->limbs[idx]; // low part
            carry = (limb_t)(res128 >> 64) + (sum < res->limbs[idx]); // hi + (0 or 1)
            // propagate the carry until it dies out
            int k = idx - 1;
            while (carry && k >= 0){
                __uint128_t temp = (__uint128_t)res->limbs[k] + carry;
                res->limbs[k] = (limb_t)temp;
                carry = temp >> 64;
                k--;
            }
        }
    }
}


static __m256i _mm256_mulhi_epi64(arbfloat_t *op1, arbfloat_t *op2, int i, int j){
    return _mm256_set_epi64x(
        (limb_t)((((__uint128_t) op1->limbs[i]) * op2->limbs[j+3]) >> 64),
        (limb_t)((((__uint128_t) op1->limbs[i]) * op2->limbs[j+2]) >> 64),
        (limb_t)((((__uint128_t) op1->limbs[i]) * op2->limbs[j+1]) >> 64),
        (limb_t)((((__uint128_t) op1->limbs[i]) * op2->limbs[j+0]) >> 64)
    ); // does op1->limbs[i] * 4 op2 limbs, and stores the high bits
}


// for (int i=(op1->length % 4) -1; i >= 0; i--){
/*
    if (op1->length < op2->length){ // make sure op1 is always longer for max loop speed
        arbfloat_t *temp = op1;
        op1 = op2;
        op2 = temp;
    }
        
*/
#ifdef __AVX512F__

void multsimd(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t *res){
    // Naive, schoolbook multiplication, accelerated with AVX2
    __m256i mask = _mm256_set_epi64x(0, -1, -1, -1);
    __m256i next_carry_mask = _mm256_set_epi64x(-1, 0, 0, 0);
    __m256i carryover = SIMD_ALL_ZEROES;
    __m256i sum, carry;
    limb_t limbsum, limbcarry;
    // ensure op2 (the accelerated one) is the longer operand
    if (op1->length > op2->length){
        arbfloat_t *temp = op1;
        op1 = op2;
        op2 = temp;
    } 
    for (int i=op1->length-1; i>=0; i--){
        limb_t op1limb = op1->limbs[i];
        __m256i op1_group = _mm256_set1_epi64x(op1limb); // limb_t is already a long long
        for (int j=op2->length-4; j>=0; j-=4){ // Changed j: since it's a 4-limb vector, needs space when starting
            int idx = i + j + 1;
            // get the high and low bits for 4 limbs at a time
            __m256i los = _mm256_mullo_epi64(op1_group, _mm256_loadu_si256((const __m256i *)(op2->limbs + j)));
            __m256i his = _mm256_mulhi_epi64(op1, op2, i, j);
            __m256i rescur = _mm256_loadu_si256((const __m256i *)(res->limbs + idx));
            // add the 4 limbs to the current results
            sum = _mm256_add_epi64(rescur, los);
            // <=> sum < rescur (compute extra carry of 1s)
            __m256i cmp_gt = _mm256_cmpgt_epu64(rescur, sum);
            cmp_gt = _mm256_srli_epi64(cmp_gt, 63);
            // add high bits to carry (hi + 1s can't overflow twice) 
            carry = _mm256_add_epi64(his, cmp_gt);
            carry = _mm256_permute4x64_epi64(carry, 0x39);
            // separate lanes 0-2 and 3
            carryover = _mm256_and_si256(carry, next_carry_mask);
            carry = _mm256_and_si256(carry, mask);
            // propagate carries until all are disappeared from carry
            while (!_mm256_testz_si256(carry, carry)){ 
                __m256i new_result = _mm256_add_epi64(sum, carry);
                cmp_gt = _mm256_cmpgt_epu64(sum, new_result);
                // Replace sum with the carried version
                sum = new_result;
                cmp_gt = _mm256_permute4x64_epi64(cmp_gt, 0x39); 
                cmp_gt = _mm256_srli_epi64(cmp_gt, 63);
                // get rid of the last limb 
                carry = _mm256_and_si256(cmp_gt, mask);
                // propagate lost limb carry
                carryover = _mm256_add_epi64(carryover, _mm256_and_si256(cmp_gt, next_carry_mask));
            }
            _mm256_storeu_si256((__m256i *)(res->limbs + idx), sum);
            // propagate anything left
            limb_t carry_temp = _mm256_extract_epi64(carryover, 3);
            int k = idx - 1;
            while (carry_temp && k >= 0){
                __uint128_t temp = (__uint128_t)res->limbs[k] + carry_temp;
                res->limbs[k] = (limb_t)temp;
                carry_temp = temp >> 64;
                k--;
            }
        }
        for (int j=(op2->length % 4) -1; j >= 0; j--){
            int idx = i + j + 1;
            // get the high and low parts of the multiplication
            __uint128_t res128 = ((__uint128_t)op1limb * (__uint128_t)op2->limbs[j]); // 81
            limbsum = (limb_t)res128 + res->limbs[idx]; // low part
            limbcarry = (limb_t)(res128 >> 64) + (limbsum < res->limbs[idx]); // hi + (0 or 1)
            // propagate the carry until it dies out
            int k = idx - 1;
            while (limbcarry && k >= 0){
                __uint128_t temp = (__uint128_t)res->limbs[k] + limbcarry;
                res->limbs[k] = (limb_t)temp;
                limbcarry = temp >> 64;
                k--;
            }
        }
    }
}

#endif

static void copy_arbfloat_with_offset(arbfloat_t* dest, const arbfloat_t *src, int offset, int limbs) {
    dest->length = (limbs) ? limbs : src->length;
    dest->capacity = (limbs) ? limbs * 2 : src->capacity;
    dest->exponent = src->exponent;
    dest->sign = src->sign;
    dest->type = src->type;

    // Allocate memory for the limbs
    dest->limbs = malloc(dest->capacity * sizeof(limb_t));

    // Copy the limbs
    memcpy(dest->limbs, src->limbs + offset, dest->length * sizeof(limb_t));
    
    return;
}

static void init_mult_limbs(arbfloat_t *dest, int l){
    dest->length = l;
    dest->capacity = l * 2;
    // Allocate memory for the limbs
    dest->limbs = malloc(dest->capacity * sizeof(limb_t));
}

static void add_limbwise_dp(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, bool *ov){
    arbfloat_t *max = (op1->length >= op2->length) ? op1 : op2;
    arbfloat_t *min = (op1->length >= op2->length) ? op2 : op1;
    int diff = max->length - min->length;
    if (!diff){ // same precision
        bool overflow = false, temp = false; // indicates if need to carry something for next add. At start there is no previous carry
        for (int i=max->length-1; i >= 0; i--){
            limb_t result = op1->limbs[i] + op2->limbs[i];
            temp = (result < op1->limbs[i]);
            result += overflow;
            overflow = (temp || result < op1->limbs[i]);
            res->limbs[i] = result;
        }
        if (ov) *ov = overflow;
        return;
    }
    else{
        bool overflow = false, temp = false; // indicates if need to carry something for next add. At start there is no previous carry
        int resoffset = diff + (res->length - max->length); // example 123 + 45, res = 0000, so diff = 1, resoffset = 1 + (4-3) = 2, so minlength goes 1, 0 (for 45's 5 and 4), max goes 2, 1, 0 (3, 2, 1) and res goes 3, 2, 1, 0
        for (int i=min->length-1; i >= 0; i--){
            limb_t op1limb = max->limbs[i+diff];
            limb_t result = op1limb + min->limbs[i];
            temp = (result < op1limb);
            result += overflow;
            overflow = (temp || result < op1limb);
            res->limbs[i+resoffset] = result;
        }
        for (int i=diff-1, newoffset=resoffset-diff; i>=0; i--){ // any remaining limbs which precede the minimum, e.g. the 1 in 123 + 45. Needs to be added to whatever carried out
            limb_t result = max->limbs[i];
            result += overflow;
            overflow = (result < max->limbs[i]);
            res->limbs[i+newoffset] = result; //e.g. if reslimbs are 0-3 big endian, and 123+45 has the 1 left over, index 0 for 123 (1) goes in index 1 for res, so i+offset=0+1
        }
        if (overflow){ // if there still exists an overflow it'll go into the '+1' limb. This limb can either be 0 or 1
            res->limbs[0]++; // overflow goes into the extra allocated limb
        }
        if (ov) *ov = false;
        return;
    }
    
}

static void sub_limbwise_dp(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res){
    arbfloat_t *max = (op1->length >= op2->length) ? op1 : op2;
    arbfloat_t *min = (op1->length >= op2->length) ? op2 : op1;
    int diff = max->length - min->length;
    if (!diff){ // same precision
        bool borrow = false, temp = false;
        for (int i=op1->length-1; i>=0; i--){
            limb_t result = op1->limbs[i] - op2->limbs[i];
            temp = (result > op1->limbs[i]);
            result -= borrow;
            borrow = (temp || result > op1->limbs[i]);
            res->limbs[i] = result;
        }
    }
    else{
        bool borrow = false, temp = false; // indicates if need to carry something for next add. At start there is no previous carry
        int resoffset = diff + (res->length - max->length); // example 123 + 45, res = 0000, so diff = 1, resoffset = 1 + (4-3) = 2, so minlength goes 1, 0 (for 45's 5 and 4), max goes 2, 1, 0 (3, 2, 1) and res goes 3, 2, 1, 0
        for (int i=min->length-1; i >= 0; i--){
            limb_t op1limb = max->limbs[i+diff];
            limb_t result = op1limb - op2->limbs[i];
            temp = (result > op1limb);
            result -= borrow;
            borrow = (temp || result > op1limb);
            res->limbs[i+resoffset] = result;
        }
        for (int i=diff-1, newoffset=resoffset-diff; i>=0; i--){ // any remaining limbs which precede the minimum, e.g. the 1 in 123 + 45. Needs to be added to whatever carried out
            limb_t result = 0ULL - max->limbs[i];
            result -= borrow;
            borrow = (result > max->limbs[i]);
            res->limbs[i+newoffset] = result; //e.g. if reslimbs are 0-3 big endian, and 123+45 has the 1 left over, index 0 for 123 (1) goes in index 1 for res, so i+offset=0+1
        }
        if (borrow){ // if there still exists an overflow it'll go into the '+1' limb. This limb can either be 0 or INT_MAX
            res->limbs[0] = 18446744073709551615ULL; // overflow goes into the extra allocated limb
        }
        return;
    }
    
}

static int cmp_limbs(arbfloat_t *x, arbfloat_t *y){
    for (int i=0; i < (x->length >= y->length) ? x->length : y->length; i++){
        limb_t xi = (i < x->length) ? x->limbs[i] : 0;
        limb_t yi = (i < y->length) ? y->limbs[i] : 0;
        if (xi > yi) return 1;
        if (xi < yi) return -1;
    }
    return 0;
}
static int max_len(arbfloat_t *x, arbfloat_t *y){
    if (x->length >= y->length) return x->length;
    return y->length;
}

static void sub_limbwise_dp_with_offset(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, int offset){
    bool borrow = false, temp = false;
    for (int i=op1->length-1-offset; i>=0; i--){
        limb_t op2limb = (i+offset < 0) ? 0 : op2->limbs[i+offset];
        limb_t result = op1->limbs[i] - op2limb;
        temp = (result > op1->limbs[i]);
        result -= borrow;
        borrow = (temp || result > op1->limbs[i]);
        res->limbs[i] = result;
    }
}

static void add_limbwise_dp_with_offset(arbfloat_t *op1, arbfloat_t *op2, arbfloat_t* res, int offset, bool *ov){
    bool overflow = false, temp = false;
    for (int i=op1->length-1-offset; i>=0; i--){
        limb_t op2limb = (i+offset < 0) ? 0 : op2->limbs[i+offset];
        limb_t result = op1->limbs[i] + op2limb;
        temp = (result < op1->limbs[i]);
        result += overflow;
        overflow = (temp || result < op1->limbs[i]);
        res->limbs[i] = result;
    }
    *ov = overflow;
}

/*
void multkaratsuba(arbfloat_t *x, arbfloat_t *y,arbfloat_t *res){
    // xy = 2^2m(z1) + 2^m(z2 - z1 - z0) + z0
    // x = (2^m)x1 + x0, y = (2^m)y1 + y0
    
    if (ARBINTR_LENGTH(x) <= KARATSUBA_THRESHOLD){
        multlimbwise(x, y, res);
        return;
    }
    // x and y are the same length, so can use length(x) unconditionally
    int m = ARBINTR_LENGTH(x) / 2; 

    arbfloat_t x1, x0, y1, y0;
    copy_arbfloat_with_offset(&x1, x, 0, m);
    copy_arbfloat_with_offset(&x0, x, m, ARBINTR_LENGTH(x) - m);
    copy_arbfloat_with_offset(&y1, y, 0, m);
    copy_arbfloat_with_offset(&y0, y, m, ARBINTR_LENGTH(y) - m);
*/
    /*if (ARBINTR_LENGTH(x) % 2 == 1){ // if length was odd, they won't be split evenly, so fix for addition
        append_limb(&x1, 0); // since m is floor(length/2), x1 and y1 will have 1 fewer limb, e.g. m = 2 when len = 5
        append_limb(&y1, 0);
    } can remove this part since i made add_limbwise_dp (different precision) */
/*
    arbfloat_t z0, z1, z2, x0x1, y0y1; // can leave type, exponent, sign unused
    init_mult_limbs(&z0, (x0.length + y0.length)); // x0 * y0 needs len(x0) + len(y0) size
    init_mult_limbs(&z1, (x1.length + y1.length)); // same for x1 * y1
    
    init_mult_limbs(&x0x1, (x0.length + 1)); // + 1 limb for overflow
    init_mult_limbs(&y0y1, (y0.length + 1)); // + 1 limb for overflow
    add_limbwise_dp(&x0, &x1, &x0x1, NULL);
    add_limbwise_dp(&y0, &y1, &y0y1, NULL);

    init_mult_limbs(&z2, (x0x1.length + y0y1.length)); // length cannot exceed length of both combined
    
    multkaratsuba(&x0, &y0, &z0); // z0 = x0 * y0
    multkaratsuba(&x1, &y1, &z1); // z1 = x1 * y1
    multkaratsuba(&x0x1, &y0y1, &z2); // z2 = (x0 + x1)(y0 + y1)

    // xy = 2^2m(z1) + 2^m(z2-z1-z0) + z0
    arbfloat_t z2z1z0, z1z0; // z2 length > z0 length >= z1 length, z2 > z0 ?= z1
    z2z1z0.sign = 0;
    init_mult_limbs(&z1z0, max_len(&z0, &z1) + 1);
    add_limbwise_dp(&z0, &z1, &z1z0, NULL);

    int cmp = cmp_limbs(&z2, &z1z0);
    init_mult_limbs(&z2z1z0, max_len(&z2, &z1z0));
    
    if (cmp == 1){ // if z2 > z1+z0, sum will have + 2^m(z2 - z1 - z0)
        sub_limbwise_dp(&z2, &z1z0, &z2z1z0);
    } else if (cmp == -1){ // if z2 <= z1+z0, sum will have -2^m(z2 - z1 - z0)
        sub_limbwise_dp(&z1z0, &z2, &z2z1z0);
        z2z1z0.sign = 1; // negative

    // res += (z1 << 2m)
    add_limbwise_with_offset(res, &z1, res, 2*m);

    // res += z0
    add_limbwise_with_offset(res, &z0, res, 0);

    // res += or -= ((z2 - z1 - z0) << m)
    if (cmp){ // if z1z2z0 is nonzero
        if (z2z1z0.sign == 1) // if z2 - z1 - z0 < 0
        {
            sub_limbwise_with_offset(res, &z2z1z0, res, m);
        } else {
            add_limbwise_with_offset(res, &z2z1z0, res, m);
        }
    }
    return;
}
*/
/*
static void append_zero_limbs(arbfloat_t *x, int n){
    bool c = false;
    while (x->length + n > x->capacity) {
        c = true;
        x->capacity *= 2;   
    }
    if (c) {
        limb_t *newlimbs = realloc(x->limbs, x->capacity * sizeof(limb_t));
        if (newlimbs == NULL) {
            fprintf(stderr, "Could not allocate memory for new limbs\n");
            exit(EXIT_FAILURE);
        }
        x->limbs = newlimbs;
    }
    x->length += n;
}
*/