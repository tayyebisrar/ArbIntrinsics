#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <assert.h>
#include <time.h>
#include "fphelpers.h"
#include "testhelpers.h"
#include "operations.h"

#define NITER 10000000 // 10 million

static void seq_right_shift_help(arbfloat_t *x, int n){
    limb_t mask = (((limb_t)9223372036854775808ULL >> (64-n-1)) - 1); // mask of k 0s then all 1s
    for (int i=0; i<ARBINTR_LENGTH(x) - 1; i++){
        x->limbs[i] = (x->limbs[i] >> n);
        x->limbs[i] |= (x->limbs[i+1] | mask);
    }
    x->limbs[ARBINTR_LENGTH(x)-1] >>= n;
}


static void seq_right_shift(arbfloat_t *x, int m, uint8_t *grt){
    if (m<=0) return;
    
    if (grt != NULL) {
        *grt = 0;
        if (m <= 64){ // likely
            if (m >= 1) *grt |= ((x->limbs[0] >> (m-1)) & 1) << 2; // G
            if (m >= 2) *grt |= ((x->limbs[0] >> (m-2)) & 1) << 1; // R
            if (m >= 3) *grt |= (x->limbs[0] & ((1ULL << (m-2)) - 1)) ? 1 : 0; // T
        }
        else{
            int limb_idx = (m / 64), bit_idx = (m % 64);
             // G
            if (bit_idx){ // bit_idx > 0, 
                *grt |= ((x->limbs[limb_idx] >> (bit_idx-1)) & 1) << 2;
            }
            else {
                *grt |= ((x->limbs[limb_idx-1] >> 63) & 1) << 2; // G is the top bit of prev limb (64th) 
            }

            // R
            if (bit_idx <= 1){
                *grt |= ((x->limbs[limb_idx-1] >> (62 + bit_idx)) & 2) << 1; // if G is at the beginning of the next (bitidx = 0, R is the 2nd along (reach by rightshifting by 62))
                // if G is at the end of the first limb (bitidx = 1), R is at the beginning of the next, reach by rightshifting 63 (62 + 1).
            }
            else{
                *grt |= ((x->limbs[limb_idx] >> (bit_idx - 2)) & 2) << 1;
                // R is in the same limb as G
            }

            // T
            // Loop through totally-lost limbs and lost bits, except G and R, to find anything nonzero
            // bit_idx = 2 is easiest - just check all preceding limbs
            // bit_idx > 2 means checking same limb without G and R
            // bit_idx < 2 means checking prev limb without G and R

            // combine = 2 and <= 2 using << 2-bit_idx (if = 2, checks whole limb (<< 0), else checks without G (=1) or GR (=0))
            if (bit_idx <= 2){
                if ((x->limbs[limb_idx - 1]) << (2-bit_idx)){
                    *grt |= 1;
                }
                if (!(*grt & 1)){ // still need to find a 1?
                    for (int i=limb_idx-2; i>=0; i--){
                        if (x->limbs[i]){
                            *grt |= 1;
                            break;
                        }
                    }
                }
            }
            else{
                *grt |= (x->limbs[limb_idx] & ((1ULL << (bit_idx-2)) - 1)) ? 1 : 0;

                if (!(*grt & 1)){
                    for (int i=limb_idx-1; i>=0; i--){
                        if (x->limbs[i]){
                            *grt |= 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    for (int i = 0, n = 32; i < m / 32; i++){
        x->exponent += n;   
        // limb_t lost_bits is a 64-bit integer, while grt is a pointer to a uint8_t, which do not match in size
        limb_t mask = (((limb_t)9223372036854775808ULL >> (64-n-1)) - 1); // mask of k 0s then all 1s
        for (int i=0; i<ARBINTR_LENGTH(x) - 1; i++){
            x->limbs[i] = (x->limbs[i] >> n);
            x->limbs[i] |= (x->limbs[i+1] | mask);
        }
        x->limbs[ARBINTR_LENGTH(x)-1] >>= n;
    }
    if ((m % 32)){
        int n = m % 32;
        x->exponent += n;
        // limb_t lost_bits is a 64-bit integer, while grt is a pointer to a uint8_t, which do not match in size
        limb_t mask = (((limb_t)9223372036854775808ULL >> (64-n-1)) - 1); // mask of k 0s then all 1s
        for (int i=0; i<ARBINTR_LENGTH(x) - 1; i++){
            x->limbs[i] = (x->limbs[i] >> n);
            x->limbs[i] |= (x->limbs[i+1] | mask);
        }
        x->limbs[ARBINTR_LENGTH(x)-1] >>= n;
    }

    // MS limb may now have garbage values at the start so need to clear them
    // x->limbs[ARBINTR_LENGTH(x)-1] &= ;
}

static void test_shifts(){
    arbfloat_t x;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    for (int i = 0; i < NITER; i++){
        generate_random_arbfloat(&x, 121);
        seq_right_shift(&x, 32, NULL);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    double elapsed = (end.tv_sec - start.tv_sec) +
                     (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Sequential shift: %d iterations took %.6f seconds\n", NITER, elapsed);
}

int main(void){
    test_shifts();
    return 0;
}