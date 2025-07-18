#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fphelpers.h>
#include <testhelpers.h>
#include <operations.h>

int main(void){
    srand(42);
    arbfloat_t x, y, res;
    generate_random_arbfloat(&x, 4);
    generate_random_arbfloat(&y, 4);  
    addnaive(&x, &y, &res);
    printf("x\n");
    display_arbfloat_limbs(&x, 1);
    printf("y\n");
    display_arbfloat_limbs(&y, 1);
    printf("x + y = \n");
    display_arbfloat_limbs(&res, 1);
}