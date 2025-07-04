#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>

void normalise(arbfloat_t *x){
    // detect first 1 (if any) - LOD detector
    // shift to it
    // return
}


void left_shift(arbfloat_t *x, int n){
    for (int i=0; i<x->length; i+=256){
        // chunks of 4 limbs
    }
    // load((x->limbs[i], x->limbs[i+1], x->limbs[i+2], x->limbs[i+3]));
    __m256i curr4 = _mm256_loadu_si256(*(((x->limbs) + 4*i)));
    // do standard shift to swap ks to the front - goes from kl to lk
    __m256i l_part = _mm256_slli_epi64(curr4, n); // l <-
    __m256i k_part = _mm256_srli_epi64(curr4, (64-n)); // k ->
    __m256i swapped = _mm256_or_si256(l_part, k_part);
    // permute the 64-bit integers across lanes, i.e. shift them around within the array such that the 0 is at the end (separate copy)
    __m256i permuted = _mm256_permute4x64_epi64(swapped, 198); // (0b11000110), ipad for logic
    // blend swapped and permuted to create shifted pattern
    // mask/separate the first 64-n bits from each integer in first (permuted) set
    // and combine with the last n bits from each integer in second (permuted and shifted) set
    __m256i lfinals = _mm256_slli_epi64(swapped, n); // l <- 
    __m256i kfinals = _mm256_srli_epi64(permuted, 64-n); // k ->
    __m256i swapped = _mm256_or_si256(lfinals, kfinals);

    // take the first k bits from the next limb (if it exists) and replace
    // the last k bits from the last limb

}


void right_shift(arbfloat_t *x, int n, uint8_t *stickyr){ // for doing (R+T)
    // similar to left_shift but backwards
    // take the last 2 bits, can be done at start, and do last OR secondlast = stickyr
} 