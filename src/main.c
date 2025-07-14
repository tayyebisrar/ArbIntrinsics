#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fphelpers.h>
#include <testhelpers.h>

int main(void){
    printf("\n");
    srand(time(NULL));
    printf("\n Before: \n");
    srand(42);
    arbfloat_t x;
    generate_random_arbfloat(&x, 4); // Generate a random arbfloat with 4 limbs
    display_arbfloat_limbs(&x, 0);
    uint8_t grt;
    printf("\n");
    for (int i=1; i<5; i++){
        printf("Right shift by %d:\n", i);
        right_shift(&x, i, &grt); // Right shift by i bits
        display_arbfloat_limbs(&x, 0);
        printf("\nGRT: ");
        for (int i=4; i>=1; i/=2){
            if (grt & i) {
                printf("1");
            } else {
                printf("0");
            }
        }
        printf("\n");
        printf("\n\n");
    }
    normalise(&x); // Normalise the arbfloat
    printf("Normalised:\n");
    display_arbfloat_limbs(&x, 0);
}
