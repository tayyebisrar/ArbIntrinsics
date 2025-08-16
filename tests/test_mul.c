// Tests the addnaive() function's effective addition

#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <assert.h>
#include "testhelpers.h"
#include "operations.h"
#include "mul.h"

static void test_multnaive(){
    // Operands are not normalised
    // Tests for correct multiplication only
    arbfloat_t x, y, res, expected;
    res.type = ARBINTR_TYPE_NORMAL; // set type to normal
    limb_t limbs_x[] = {
        0x1234567890abcdef,
        0xfedcba0987654321
    }; // 00010010001101000101011001111000100100001010101111001101111011111111111011011100101110100000100110000111011001010100001100100001
    new_arbfloat(&x, limbs_x, 2, 1, 0, ARBINTR_TYPE_NORMAL); // x, limbs, length, exponent, sign, normal

    limb_t limbs_y[] = {
        0x1111111111111111,
        0x2222222222222222
    };// 00010001000100010001000100010001000100010001000100010001000100010010001000100010001000100010001000100010001000100010001000100010
    new_arbfloat(&y, limbs_y, 2, 1, 0, ARBINTR_TYPE_NORMAL); // y, limbs, length, exponent, sign, normal

    
    limb_t expected_limbs[] = {
        0x0136b06e700b7421,
        0x1234566f12abcdf0
        //0x0e904523c289718c,
        //0xde04b3febabf4c62
    }; // 0x0136b06e700b74211234566f12abcdf00e904523c289718cde04b3febabf4c62 with the last half gone
    new_arbfloat(&expected, expected_limbs, 2, 2, 0, ARBINTR_TYPE_NORMAL);
    multnaive(&x, &y, &res, RNDZ);
    assert((arbfloat_equal(&res, &expected) == 1) && "Multnaive failed");
    printf("Test multnaive passed.\n");
}

int main(void){
    test_multnaive();
    return 0;
}