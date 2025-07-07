#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <string.h>

void normalise(arbfloat_t *x){
    // detect first 1 (if any) - LOD detector
    int NElimb = 0; // index for first non-zero limb
    for (int i=0; i<x->length; i++){
        if (x->limbs[i] != 0){
            NElimb = i;
            break;
        } 
    }
    if (!NElimb && x->limbs[0] == 0){
        // number is all zeroes
        return;
    }
    for (int i=0; i<NElimb; i++){ // for each limb 
        left_shift(&x, 63);
        left_shift(&x, 1);
        x->length--;
    }

    (x->limbs[NElimb]);
    // shift to it
    // return
}


void left_shift(arbfloat_t *x, int n){
    // CAN BE HERE TO SAVE COPYING -- limb_t mask = ~(((limb_t)1 << n) - 1); // bitmask, starts as all 1s and k 0s at the end
    // Runs from [0, length-4) to stop lengths divisible by 4 from running the last block SIMD
    for (int i=0; i<x->length - 4; i+=4){ // process 4 limbs (256 bits) at a time
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
        _mm256_storeu_si256((__m256i *)(x->limbs + i), finalshifted);

        // take the first k bits from the next limb (if it exists) and replace
        // the last k bits from the last limb, i.e. for (int i=63-k+1; i<=63; i++){limb4[i] = limb5[i]}
        limb_t mask = ~(((limb_t)1 << n) - 1); // bitmask, starts as all 1s and k 0s at the end
        x->limbs[i+3] &= mask; // clear the last k bits
        limb_t limb5 = x->limbs[i+4] >> (64-n); // next limb shifted down, going from kklll to 000..kk
        x->limbs[i+3] |= limb5; // new bitmask to set limb4's last k bits via bitwise AND. i.e., 111..0 | 000...k = 111...k
        
    }
    // process any remaining limbs that are not a multiple of 4
    limb_t mask = ~(((limb_t)1 << n) - 1);
    int seqiter = (x->length / 4); // 1
    if (x->length % 4 == 0){
        seqiter--;
    }
    for (int i=seqiter*4; i<x->length; i++){    
        x->limbs[i] <<= n;
        if (i != x->length - 1){ // if k bits need to be shifted, from the next one in
            x->limbs[i] &= mask;
            limb_t limb2 = x->limbs[i+1] >> (64-n);
            x->limbs[i] |= limb2;
        }
    }

    // remove any now-empty limb(s) from the end of the array
    while (x->length > 0 && x->limbs[x->length - 1] == 0) {
        //free((limb_t *)(x->limbs[x->capacity - 1])); wont work because it's a contiguous block of memory
        x->length--; // logically get rid of it like a stack pointer
    }
}


void right_shift(arbfloat_t *x, int n, uint8_t *stickyr){ // for doing (R+T)
    // take the last 2 bits, can be done at start, and do last OR secondlast = stickyr
    // Do the first 4 separately, since otherwise it will access limbs[-1]]

    limb_t mask = (((limb_t)9223372036854775808 >> (n-1)) - 1);
    for (int i=1, a=min(4, x->length); i<a ; i++){
        x->limbs[i] >>= n;
        x->limbs[i] &= mask; // clear the first k bits
        if (i != 0){
            limb_t limb5 = x->limbs[i-1] << (64-n); // next limb shifted up, going from lllkk to kk000..00
            x->limbs[i] |= limb5; // new bitmask to set limb4's last k bits via bitwise AND. i.e., 000..1010 | kk...00 = kk...1010
        }
    }
    for (int i=4; i<x->length; i+=4){
        __m256i curr4 = _mm256_loadu_si256((const __m256i *)(x->limbs + i));
        // do standard shift to swap ks to the back - goes from lk to kl
        __m256i l_part = _mm256_srli_epi64(curr4, n); // l ->
        __m256i k_part = _mm256_slli_epi64(curr4, (64-n)); // k <-
        __m256i swapped = _mm256_or_si256(l_part, k_part);
        __m256i permuted = _mm256_permute4x64_epi64(swapped, 108); // (0b01101100), ipad for logic
        // mask/separate the last 64-n bits from each integer in first (permuted) set
        // and combine with the first n bits from each integer in second (permuted and shifted) set
        __m256i lfinals = _mm256_srli_epi64(swapped, n); // l ->
        __m256i kfinals = _mm256_slli_epi64(permuted, 64-n); // k <-
        __m256i finalshifted = _mm256_or_si256(lfinals, kfinals);
        _mm256_storeu_si256((__m256i *)(x->limbs + i), finalshifted);

        // take the last k bits from the previous limb (if it exists) and replace
        // the first k bits from the first limb, i.e. for (int i=63-k+1; i<=63; i++){limb4[i] = limb5[i]}
        limb_t mask = (((limb_t)9223372036854775808 >> (n-1)) - 1); // bitmask, starts with 0 ks and then all 1s
        x->limbs[i] &= mask; // clear the first k bits
        limb_t limb5 = x->limbs[i-1] << (64-n); // prev limb shifted up, going from lllkk to kk000..00
        x->limbs[i] |= limb5; // new bitmask to set limb4's last k bits via bitwise AND. i.e., 000..1010 | kk...00 = kk...1010
    }
    for (int i=(x->length/4)*4; i<x->length; i++){    
        x->limbs[i] >>= n;
        x->limbs[i] &= mask; // clear the first k bits
        limb_t limb5 = x->limbs[i-1] << (64-n); // next limb shifted up, going from lllkk to kk000..00
        x->limbs[i] |= limb5; // new bitmask to set limb4's last k bits via bitwise AND. i.e., 000..1010 | kk...00 = kk...1010
    }

    // first limb will have garbage values at the start so need to clear them
    x->limbs[0] &= mask;
}


