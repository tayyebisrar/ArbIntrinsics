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
            tobinary(a->limbs[i]);
            printf("\nvs\n");
            tobinary(b->limbs[i]);
            printf("\n");
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
    assert((arbfloat_equal(&res, &expected) == 1) && "Regular addition failed");
    printf("Test regular addition passed.\n");
}

static void testzero(){
    arbfloat_t x, y, res, expected;
    limb_t limbs_x[] = {0x0}; // 0
    new_arbfloat(&x, limbs_x, 1, 0, 0, ARBINTR_TYPE_ZERO); // x, limbs, length, exponent, sign, zero

    limb_t limbs_y[] = {0x0}; // 0
    new_arbfloat(&y, limbs_y, 1, 0, 0, ARBINTR_TYPE_ZERO); // y, limbs, length, exponent, sign, zero

    arbintr_add(&x, &y, &res, RNDN);
    new_arbfloat(&expected, limbs_x, 1, 0, 0, ARBINTR_TYPE_ZERO); // expected result is also zero
    assert((res.type == ARBINTR_TYPE_ZERO) && "Zero + zero addition failed");
    printf("Zero + zero addition passed.\n");
}

static void testzeroadds(){
    arbfloat_t x, y, res, expected;
    limb_t limbs_x[] = {0x0}; // 0
    new_arbfloat(&x, limbs_x, 1, 0, 0, ARBINTR_TYPE_ZERO); // x, limbs, length, exponent, sign, zero

    limb_t limbs_y[] = {0x1234567890ABCDEF}; // 1 limb
    new_arbfloat(&y, limbs_y, 1, 0, 0, ARBINTR_TYPE_NORMAL); // y, limbs, length, exponent, sign, normal

    limb_t expected_limbs[] = {0x1234567890ABCDEF};
    new_arbfloat(&expected, expected_limbs, 1, 0, 0, ARBINTR_TYPE_NORMAL); // expected result
    
    arbintr_add(&x, &y, &res, RNDN);

    bool first_test = (arbfloat_equal(&res, &expected) == 1);
    
    arbintr_add(&y, &x, &res, RNDN);

    bool second_test = (arbfloat_equal(&res, &expected) == 1);

    assert((first_test) && "Zero + y addition failed");
    printf("Zero + y addition passed.\n");
    assert((second_test) && "x + zero addition failed");
    printf("x + zero addition passed.\n");
}

static void testinf(){
    arbfloat_t x, y, res, expected;
    limb_t limbs_x[] = {0x0}; // 0
    new_arbfloat(&x, limbs_x, 1, 0, 0, ARBINTR_TYPE_INF); // x, limbs, length, exponent, sign, infinity

    limb_t limbs_y[] = {0x0}; // 0
    new_arbfloat(&y, limbs_y, 1, 0, 0, ARBINTR_TYPE_NORMAL); // y, limbs, length, exponent, sign, normal

    arbintr_add(&x, &y, &res, RNDN);

    bool first_test = (res.type == ARBINTR_TYPE_INF);
    
    x.sign = 1; // x = -inf

    arbintr_add(&res, &x, &res, RNDN); // Inf - Inf = NaN

    bool second_test = (res.type == ARBINTR_TYPE_NAN);

    assert((first_test) && "Normal + Inf addition failed");
    printf("Normal + Inf addition passed.\n");
    assert((second_test) && "Inf-Inf (NaN) addition failed");
    printf("Inf-Inf (NaN) addition passed.\n");
}

static void testround(){
    // Same values as test_regular but with increase exponent difference for rounding
    arbfloat_t x, y, res, expected;
    limb_t limbs_x[] = {
        0x3ed3756845cd0190,
        0x3dc11ddc1f50324f,
        0x349f74fe7e9906e9,
        0x52c45ab909201bbb
    }; // srand(42) yields this
    new_arbfloat(&x, limbs_x, 4, 8, 0, ARBINTR_TYPE_NORMAL); // x, limbs, length, exponent, sign, normal

    limb_t limbs_y[] = {
        0x686a22cb4cdc380f,
        0x1b8023e6640773f6,
        0x0e830d7517cb18c0,
        0x442570d80e0d75e6
    }; 
    new_arbfloat(&y, limbs_y, 4, 0, 0, ARBINTR_TYPE_NORMAL); // y, limbs, length, exponent, sign, normal

    arbintr_add(&x, &y, &res, RNDN);
    
    limb_t expected_limbs[] = {
        0xfcef7e2c44677721,
        0x3372780016d0e70c,
        0xaab7e02fcec34808,
        0x4c2200a784b8a4c4
    }; // 0xfcef7e2c446777213372780016d0e70caab7e02fcec348084c2200a784b8a4c4
    new_arbfloat(&expected, expected_limbs, 4, 6, 0, ARBINTR_TYPE_NORMAL);
    assert((arbfloat_equal(&res, &expected) == 1) && "Rounded addition failed");
    printf("Test rounded addition passed.\n");
}

static void test_subtract(){
    // x - y
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
    new_arbfloat(&y, limbs_y, 4, 0, 1, ARBINTR_TYPE_NORMAL); // y, limbs, length, exponent, sign, normal

    arbintr_add(&x, &y, &res, RNDN);
    
    limb_t expected_limbs[] = {
        0x92e3b2d5ca57ce31,
        0xdb84538a19395546,
        0xc3fac684e29902e5,
        0x06ebfa0c1672f906
    }; // 0x92e3b2d5ca57ce31db84538a19395546c3fac684e29902e506ebfa0c1672f906
    new_arbfloat(&expected, expected_limbs, 4, 0, 0, ARBINTR_TYPE_NORMAL);
    assert((arbfloat_equal(&res, &expected) == 1) && "Regular subtraction failed");
    printf("Test regular subtraction passed.\n");
}

static void test_subround(){
    // x - y with a higher exp. difference to see rounding
    arbfloat_t x, y, res, expected;
    limb_t limbs_x[] = {
        0x3ed3756845cd0190,
        0x3dc11ddc1f50324f,
        0x349f74fe7e9906e9,
        0x52c45ab909201bbb
    }; // srand(42) yields this
    new_arbfloat(&x, limbs_x, 4, 8, 0, ARBINTR_TYPE_NORMAL); // x, limbs, length, exponent, sign, normal

    limb_t limbs_y[] = {
        0x686a22cb4cdc380f,
        0x1b8023e6640773f6,
        0x0e830d7517cb1c0,
        0x442570d80e0d75e6
    }; 
    new_arbfloat(&y, limbs_y, 4, 0, 1, ARBINTR_TYPE_NORMAL); // y, limbs, length, exponent, sign, normal

    arbintr_add(&x, &y, &res, RNDA);
    
    limb_t expected_limbs[] = {
        0xf9ac2d15ea009560,
        0xba9676e0e3b0ab6c,
        0xfa43c7c42604ef42,
        0x4a00d520c4483915
    }; // 0xf9ac2d15ea009560ba9676e0e3b0ab6cfa43c7c42604ef424a00d520c4483915
    new_arbfloat(&expected, expected_limbs, 4, 6, 0, ARBINTR_TYPE_NORMAL);
    assert((arbfloat_equal(&res, &expected) == 1) && "Rounded subtraction failed");
    printf("Test rounded subtraction passed.\n");
}




int main(void){
    test_addnaive();
    test_regular();
    testround();
    testzero();
    testzeroadds();
    testinf();
    test_subtract();
    test_subround();
    return 0;
}