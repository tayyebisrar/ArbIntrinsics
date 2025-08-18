// Tests the arbintr_mul() function's multiplcation

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

static void test_regular(){
    arbfloat_t x, y, res, expected;
    limb_t limbs_x[] = {
        0x3ed3756845cd0190,
        0x3dc11ddc1f50324f,
        0x349f74fe7e9906e9,
        0x52c45ab909201bbb
    }; // srand(42) yields this
    new_arbfloat(&x, limbs_x, 4, 2, 0, ARBINTR_TYPE_NORMAL); // x, limbs, length, exponent, sign, normal

    limb_t limbs_y[] = {
        0x686a22cb4cdc380f,
        0x1b8023e6640773f6,
        0x0e830d7517cb18c0,
        0x442570d80e0d75e6
    };
    new_arbfloat(&y, limbs_y, 4, 0, 0, ARBINTR_TYPE_NORMAL); // y, limbs, length, exponent, sign, normal

    arbintr_mul(&x, &y, &res, RNDN);
    
    limb_t expected_limbs[] = {
        0xccff9e579d57a67d,
        0x1de5a4ed45899113,
        0x534acb2da4296f93,
        0x324aaccdb574531b,
    }; // 0xccff9e579d57a67d1de5a4ed45899113534acb2da4296f93324aaccdb574531b2f8b9cc34f5e7fd62a5f0be35f4878702c2772137c23b612082c0abc10230810
    // with last half gone: 0xccff9e579d57a67d1de5a4ed45899113534acb2da4296f93324aaccdb574531b
    // G=0, T=1
    // final: 0xccff9e579d57a67d1de5a4ed45899113534acb2da4296f93324aaccdb574531b with RNDN
    new_arbfloat(&expected, expected_limbs, 4, -1, 0, ARBINTR_TYPE_NORMAL);
    assert((arbfloat_equal(&res, &expected) == 1) && "Regular multiplication failed");
    printf("Test regular multiplication passed.\n");
}

int main(void){
    test_multnaive();
    test_regular();
    return 0;
}