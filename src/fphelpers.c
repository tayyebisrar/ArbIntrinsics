#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <immintrin.h>
#include <string.h>
#include <testhelpers.h>

void left_shift(arbfloat_t *x, int n);
void right_shift(arbfloat_t *x, int n, uint8_t *grt);
void normalise(arbfloat_t *x);


void normalise(arbfloat_t *x){
    // left only
    // detect first 1 (if any) - LOD detector
    int NElimb = -1; // index of first non-zero limb
    for (size_t i = 0; i < x->length; i++){
        if (x->limbs[i] != 0){
            NElimb = i;
            break;
        } 
    }
    if (NElimb == -1) {
        // All zero
        x->length = 0;
        x->exponent = 0;
        x->sign = 0;
        return;
    }
    int leadingzeroes = __builtin_clzll(x->limbs[NElimb]);
    for (int i=0; i<NElimb; i++){
        left_shift(x, 32);
        left_shift(x, 32);
        x->exponent -= 64;
    }
    left_shift(x, leadingzeroes);
    return;
    // return
}


void left_shift(arbfloat_t *x, int n){
    x->exponent -= n; // shift exponent down
    limb_t mask = ~(((limb_t)1 << n) - 1); // bitmask, starts as all 1s and k 0s at the end
    // Runs from [0, length-4) to stop lengths divisible by 4 from running the last block SIMD
    int seqiter = (x->length / 4);
    if (x->length % 4 == 0){
        seqiter--;
    }
    for (int i=0; i<seqiter*4; i+=4){ // process 4 limbs (256 bits) at a time
        __m256i curr4 = _mm256_loadu_si256((const __m256i *)(x->limbs + i)); // load the next 256 bits of integer data
        // i.e. // load((x->limbs[i], x->limbs[i+1], x->limbs[i+2], x->limbs[i+3]));
        // do standard shift to swap ks to the front - goes from kl to lk
        __m256i l_part = _mm256_slli_epi64(curr4, n); // l <-
        __m256i k_part = _mm256_srli_epi64(curr4, (64-n)); // k ->
        // permute the 64-bit integers across lanes, i.e. shift them around within the array such that the 0 is at the end (separate copy)
        __m256i permuted = _mm256_permute4x64_epi64(k_part, 0x39); // (0b00111001), was 0b10010011 0b11000110 was backwards 
        // blend swapped and permuted to create shifted pattern
        // mask/separate the first 64-n bits from each integer in first (permuted) set
        // and combine with the last n bits from each integer in second (permuted and shifted) set
        __m256i finalshifted = _mm256_or_si256(l_part, permuted);
        // store back the 4 limbs
        _mm256_storeu_si256((__m256i *)(x->limbs + i), finalshifted);

        // take the first k bits from the next limb (if it exists) and replace
        // the last k bits from the last limb, i.e. for (int i=63-k+1; i<=63; i++){limb4[i] = limb5[i]}
        x->limbs[i+3] &= mask; // clear the last k bits
        limb_t limb5 = x->limbs[i+4] >> (64-n); // next limb shifted down, going from kklll to 000..kk
        x->limbs[i+3] |= limb5; // new bitmask to set limb4's last k bits via bitwise AND. i.e., 111..0 | 000...k = 111...k
    }
    // process any remaining limbs that are not a multiple of 4
    for (int i=seqiter*4; i<x->length; i++){
        x->limbs[i] <<= n;
        if (i != x->length - 1){ // if k bits need to be shifted, from the next one in
            x->limbs[i] &= mask;
            limb_t limb2 = x->limbs[i+1] >> (64-n);
            x->limbs[i] |= limb2;
        }
    }

    /* remove any now-empty limb(s) from the end of the array
    while (x->length > 0 && x->limbs[x->length - 1] == 0) {
        //free((limb_t *)(x->limbs[x->capacity - 1])); wont work because it's a contiguous block of memory
        x->length--; // logically get rid of it like a stack pointer
    } */
}


void right_shift(arbfloat_t *x, int n, uint8_t *grt){ // for doing (R+T)
    x->exponent += n;
    // limb_t lost_bits is a 64-bit integer, while grt is a pointer to a uint8_t, which do not match in size
    if (grt != NULL) {
        *grt = 0;
        limb_t lost_bits = x->limbs[x->length-1] & ((1ULL << n) - 1); // bits shifted out (= n <= 63 bits)
        if (n >= 1) *grt |= ((lost_bits >> (n-1)) & 1) << 2; // G
        if (n >= 2) *grt |= ((lost_bits >> (n-2)) & 1) << 1; // R
        if (n >= 3) *grt |= (lost_bits & ((1ULL << (n-2)) - 1)) ? 1 : 0; // T
    }
    // Do the first 4 separately, since otherwise it will access limbs[-1]]
    limb_t mask = (((limb_t)9223372036854775808ULL >> (n-1)) - 1); // mask of k 0s then all 1s
    
    for (int i=x->length-4; i>3; i-=4){
        __m256i curr4 = _mm256_loadu_si256((const __m256i *)(x->limbs + i)); // load the next 256 bits of integer data
        __m256i l_part = _mm256_srli_epi64(curr4, n); // l ->
        __m256i k_part = _mm256_slli_epi64(curr4, (64-n)); // k <-
        __m256i permuted = _mm256_permute4x64_epi64(k_part, 0x93);
        __m256i finalshifted = _mm256_or_si256(l_part, permuted);
        _mm256_storeu_si256((__m256i *)(x->limbs + i), finalshifted);

        x->limbs[i] &= mask; // clear the first k bits
        limb_t limb5 = x->limbs[i-1] << (64-n); // prev limb shifted up, going from lllkk to kk000..00
        x->limbs[i] |= limb5; // new bitmask to set limb4's last k bits via bitwise AND. i.e., 000..1010 | kk...00 = kk...1010
    }

    // Handle remaining limbs
    for (int i=(x->length < 4) ? x->length-1 : 3 + x->length % 4; i >= 0; i--){
        x->limbs[i] >>= n;
        if (i != 0){
            x->limbs[i] &= mask; // clear the first k bits
            limb_t limb5 = x->limbs[i-1] << (64-n);
            x->limbs[i] |= limb5;
        }
    }

    // first limb may have garbage values at the start so need to clear them
    x->limbs[0] &= mask;
}

arbfloat_t *copy_arbfloat(const arbfloat_t *src) {
    arbfloat_t *dest = malloc(sizeof(arbfloat_t));
    dest->length = src->length;
    dest->capacity = src->capacity;
    dest->exponent = src->exponent;
    dest->sign = src->sign;

    // Allocate memory for the limbs
    dest->limbs = malloc(dest->capacity * sizeof(limb_t));

    // Copy the limbs
    memcpy(dest->limbs, src->limbs, dest->length * sizeof(limb_t));
    
    return dest;
}

void append_limb(arbfloat_t *x, limb_t limb){
    // Increase length of x's limb array if needed, then append limb_t with value limb
    if (x->length >= x->capacity) {
        x->capacity *= 2;
        limb_t *newlimbs = realloc(x->limbs, x->capacity * sizeof(limb_t));
        if (newlimbs == NULL) {
            fprintf(stderr, "Could not allocate memory for new limbs\n");
            exit(EXIT_FAILURE);
        }
        x->limbs = newlimbs;
    }
    x->limbs[x->length] = limb;
    x->length++;
}