#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fphelpers.h>
#include <testhelpers.h>

int main(void){
    printf("\n Before: \n");
    srand(42);
    arbfloat_t x;
    generate_random_arbfloat(&x, 5); // Generate a random arbfloat with 4 limbs
    display_arbfloat_limbs(&x, 0);
    printf("\n");
    for (int i=1; i<5; i++){
        printf("Left shift by %d:\n", i);
        left_shift(&x, i); // Right shift by i bits
        display_arbfloat_limbs(&x, 0);
        printf("\n");
    }
}