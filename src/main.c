#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fphelpers.h>
#include <testhelpers.h>

int main(void){
    srand(time(NULL));
    arbfloat_t x;
    generate_random_arbfloat(&x, 4); // Generate a random arbfloat with 4 limbs
    display_arbfloat(&x);
}