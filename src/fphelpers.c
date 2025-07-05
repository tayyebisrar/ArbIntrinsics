#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <string.h>

void normalise(arbfloat_t *x){
    // detect first 1 (if any) - LOD detector
    // shift to it
    // return
}


void left_shift(arbfloat_t *x, int n){
    for (int i=0; i<x->length ; i+=4){ // process 4 limbs (256 bits) at a time
        __m256i curr4 = _mm256_loadu_si256((const __m256i *)(x->limbs + i)); // load the next 256 bits of integer data
        // i.e. // load((x->limbs[i], x->limbs[i+1], x->limbs[i+2], x->limbs[i+3]));
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
        __m256i finalshifted = _mm256_or_si256(lfinals, kfinals);
        // store back the 4 limbs
        _mm256_storeu_si256((const __m256i *)(x->limbs + i), finalshifted);

        // take the first k bits from the next limb (if it exists) and replace
        // the last k bits from the last limb, i.e. for (int i=63-k+1; i<=63; i++){limb4[i] = limb5[i]}
        limb_t mask = ~(((limb_t)1 << n) - 1); // bitmask, starts as all 1s and k 0s at the end
        x->limbs[i+3] &= mask; // clear the last k bits
        limb_t limb5 = x->limbs[i+4] >> (64-n); // next limb shifted down, going from kklll to 000..kk
        x->limbs[i+3] |= limb5; // new bitmask to set limb4's last k bits via bitwise AND. i.e., 111..0 | 000...k = 111...k
        
    }
    
    for (int i=(x->length / 4)*4; i<x->length; i++){    
        // process any remaining limbs that are not a multiple of 4
    }

    // remove any now-empty limb(s) from the end of the array
}


void right_shift(arbfloat_t *x, int n, uint8_t *stickyr){ // for doing (R+T)
    // similar to left_shift but backwards
    // take the last 2 bits, can be done at start, and do last OR secondlast = stickyr
} 

/*
limb_t mask = ~(((limb_t)1 << n) - 1); // bitmask, starts as all 1s and k 0s at the end
        limb_t limb5 = x->limbs[i+5] >> 64-n; // next limb shifted down, going from kklll to 000..kk
        mask |= limb5; // new bitmask to set limb4's last k bits via bitwise AND. i.e., 111..0 | 000...k = 111...k
        x->limbs[i+4] &= mask;
        */