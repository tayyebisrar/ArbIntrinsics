#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fphelpers.h>
#include <testhelpers.h>
#include <operations.h>

int main(void){
    srand(42);
    arbfloat_t x, y, res;
    uint8_t grt;
    generate_random_arbfloat(&x, 4);
    generate_random_arbfloat(&y, 4);  
    x.exponent = 2;
    y.exponent = 0;
    x.sign = 0;
    y.sign = 0;
    printf("x\n");
    display_arbfloat_limbs(&x, 1);
    printf("y\n");
    display_arbfloat_limbs(&y, 1);
    printf("x + y = \n");
    addnaive(&x, &y, &res, &grt);
    display_arbfloat_limbs(&res, 1);
    printf("GRT: %d\n", grt);
}