#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fphelpers.h>
#include <testhelpers.h>

int main(void){
    srand(time(NULL));
    arbfloat_t x;
    generate_random_arbfloat(&x, 1); // Generate a random arbfloat with 4 limbs
    display_arbfloat(&x, 1);
    left_shift(&x, 1);
    printf("Shifted left once: \n");
    display_arbfloat(&x, 1);
    left_shift(&x, 1);
    printf("Shifted left twice: \n");
    display_arbfloat(&x, 1);
}