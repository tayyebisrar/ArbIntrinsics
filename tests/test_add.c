// Tests the addnaive() function's effective addition

#include <stdio.h>
#include <stdlib.h>
#include <arbintrinsics.h>
#include <assert.h>
#include "testhelpers.h"
#include "operations.h"
#include "add.h"

static int arbfloat_equal(arbfloat_t *a, arbfloat_t *b) {
    if (a->length != b->length || a->exponent != b->exponent || a->sign != b->sign || a->type != b->type) {
        printf("Lengths: %zu vs %zu, Exponents: %d vs %d, Signs: %d vs %d, Types: %d vs %d\n",
               a->length, b->length, a->exponent, b->exponent, a->sign, b->sign, a->type, b->type);
        return 0; // Not equal
    }
    for (size_t i = 0; i < a->length; i++) {
        if (a->limbs[i] != b->limbs[i]) {
            return 0; // Not equal
        }
    }
    return 1; // Equal
}
static void test_addnaive(){
    // Operands are not normalised
    // Tests for correct addition only
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
        0x23456789a1bcdf01,
        0x20fedc2ba9876543
    };// 0010001101000101011001111000100110100001101111001101111100000001001000001111111011011100001010111010100110000111011001010100001100
    new_arbfloat(&expected, expected_limbs, 2, 1, 0, ARBINTR_TYPE_NORMAL);
    addnaive(&x, &y, &res, RNDN);
    assert((arbfloat_equal(&res, &expected) == 1) && "Addnaive failed");
    printf("Test addnaive passed.\n");
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

    arbintr_add(&x, &y, &res, RNDN);
    
    limb_t expected_limbs[] = {
        0xb1dbfc3632081f28,
        0x09424dab70a41e99,
        0x708070b789179a32,
        0xc79b6dde1946f269
    }; // 0xb1dbfc3632081f2809424dab70a41e99708070b789179a32c79b6dde1946f269
    new_arbfloat(&expected, expected_limbs, 4, 1, 0, ARBINTR_TYPE_NORMAL);
    assert((arbfloat_equal(&res, &expected) == 1) && "Regular addition 2 failed");
    printf("Test regular addition 2 passed.\n");
}

static void testzero(){
    arbfloat_t x, y, res, expected;
    limb_t limbs_x[] = {0x0}; // 0
    new_arbfloat(&x, limbs_x, 1, 0, 0, ARBINTR_TYPE_ZERO); // x, limbs, length, exponent, sign, zero

    limb_t limbs_y[] = {0x0}; // 0
    new_arbfloat(&y, limbs_y, 1, 0, 0, ARBINTR_TYPE_ZERO); // y, limbs, length, exponent, sign, zero

    arbintr_add(&x, &y, &res, RNDN);
    new_arbfloat(&expected, limbs_x, 1, 0, 0, ARBINTR_TYPE_ZERO); // expected result is also zero
    assert((res.type == ARBINTR_TYPE_ZERO) && "Zero addition failed");
    printf("Test zero addition passed.\n");
}

int main(void){
    test_addnaive();
    test_regular();
    testzero();
    return 0;
}