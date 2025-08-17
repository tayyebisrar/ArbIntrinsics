#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fphelpers.h>
#include <testhelpers.h>
#include <operations.h>
#include "add.h"

#include <time.h>


static void test_add(int limbs){
    arbfloat_t x, y, res;
    
    generate_random_arbfloat(&x, limbs);
    generate_random_arbfloat(&y, limbs);

    clock_t start_time = clock();

    for (int i = 0; i < 10000; i++) { // 10,000 iterations
        arbintr_add(&x, &y, &res, RNDZ);
    }

    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time)) / (double) CLOCKS_PER_SEC;
    printf("Time taken: %f seconds\n", cpu_time_used);
}

int main(void){
    srand(42);
    test_add(140);
}