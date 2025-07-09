#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fphelpers.h>
#include <testhelpers.h>

int main(void){
    srand(time(NULL));
    arbfloat_t x;
    generate_random_arbfloat(&x, 5); // Generate a random arbfloat with 4 limbs
    display_arbfloat_limbs(&x, 1);
    for (int i=1; i<5; i++){
        printf("Right shift by %d:\n", i);
        right_shift(&x, i, NULL); // Right shift by i bits
        display_arbfloat_limbs(&x, 1);
        printf("\n");
    }
    normalise(&x); // Normalise the arbfloat
    printf("Normalised:\n");
    display_arbfloat_limbs(&x, 1);
}