#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arbintrinsics.h>

#define MAX_EXPONENT 256
#define MIN_EXPONENT -256

void display_arbfloat_limbs(arbfloat_t *x, int oneline);
void tobinary(limb_t n);
limb_t rand64(void);
void generate_random_arbfloat(arbfloat_t *x, size_t length);
 
void display_arbfloat_limbs(arbfloat_t *x, int oneline){
    if (oneline){
        printf("Number: ");
        for (int i=0; i<x->length; i++){
            //printf("Limb %d: %lld\n", i, (x->limbs[i]));
            tobinary((x->limbs[i]));
        }
        printf("\n");
    }
    else{
        for (int i=0; i<x->length; i++){
            //printf("Limb %d: %lld\n", i, (x->limbs[i]));
            printf("Limb %d:", i);
            tobinary((x->limbs[i]));
            printf(" (%llu)", x->limbs[i]);
            printf("\n");
        }
    }
    printf("Exponent: %d\n", x->exponent);
}

void generate_random_arbfloat(arbfloat_t *x, size_t l) {
    x->length = l;
    x->capacity = l;
    x->limbs = (limb_t *)malloc(l * sizeof(limb_t));
    x->exponent = rand() % (MAX_EXPONENT - MIN_EXPONENT + 1); // =32767 % 1000 - (-1000) + 1 = 2001;
    x->sign = 0;

    for (int i=0; i<l; i++){
        x->limbs[i] = rand64();
    }
}

void tobinary(limb_t n) {
    char binary[65]; // 64 bits and null
    binary[64] = '\0';
    for (int i = 63; i >= 0; i--) {
        binary[i] = (n & 1) ? '1' : '0'; // Check the least significant bit
        n >>= 1; // Right shift to process the next bit
    }
    printf("%s", binary); // Print the binary representation
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

void display_normalised(arbfloat_t *x){
    printf("Number: ");
    if (x->sign){
        printf("-");
    }
    char binary[64]; // 63 bits and null
    binary[63] = '\0';
    limb_t *n = malloc(sizeof(x->limbs[0]));
    *n = x->limbs[0];
    for (int i = 62; i >= 0; i--) {
        binary[i] = (*n & 1) ? '1' : '0'; // Check the least significant bit
        *n >>= 1; // Right shift to process the next bit
    }
    free(n);
    printf("1.%s", binary);
    if (x->length > 1){
        printf("...");
    }
    printf(" x2^%d\n", x->exponent);
}