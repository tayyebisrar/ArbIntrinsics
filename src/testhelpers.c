#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arbintrinsics.h>

void display_arbfloat(arbfloat_t *x);
void tobinary(limb_t *n);
limb_t rand64(void);
void generate_random_arbfloat(arbfloat_t *x, size_t length);

void display_arbfloat(arbfloat_t *x){
    for (int i=0; i<x->length; i++){
        //printf("Limb %d: %lld\n", i, (x->limbs[i]));
        printf("Limb %d:", i);
        tobinary(&(x->limbs[i]));
        printf("\n");
    }
}

void generate_random_arbfloat(arbfloat_t *x, size_t length) {
    x->length = length;
    x->capacity = length;
    x->limbs = (limb_t *)malloc(length * sizeof(limb_t));
    x->exponent = 0;
    x->sign = 0;

    for (int i=0; i<length; i++){
        x->limbs[i] = rand64();
    }
}

void tobinary(limb_t *n) {
    char binary[65]; // 64 bits and null
    binary[64] = '\0';
    for (int i = 63; i >= 0; i--) {
        binary[i] = (*n & 1) ? '1' : '0'; // Check the least significant bit
        *n >>= 1; // Right shift to process the next bit
    }
    printf("%s\n", binary); // Print the binary representation
}

limb_t rand64(void) {
    // RAND_MAX is 32767, using 16 bits
    // To generate a 64-bit random number we can combine multiple calls to rand()
    limb_t result = 0;
    for (int i = 0; i < 4; i++) { // 4 calls to rand() to fill 64 bits
        result |= ((limb_t)rand() << (i * 16)); // Shift each 16-bit chunk into the correct position
    }
    return result;
}