# ArbIntrinsics

ArbIntrinsics is a high-precision arbitrary-precision floating-point arithmetic library implemented in C, with SIMD (AVX2/AVX512) acceleration for limb-wise operations. This library allows operations such as addition, subtraction, and multiplication on arbitrary-precision floating-point numbers using both naive and SIMD-optimised methods.

It was initially developed as a first-year internship project and demonstrates advanced techniques in numerical computing, low-level optimisation, and vectorised programming.
It was tested on DiRAC COSMA's supercomputing facility.

Actively under a complete refactor to little-endian, so at the current moment tests do not work.

Current development is unfortunately sporadic due to other commitments.

## Requirements

- C compiler (GCC/Clang or MSVC)
- CPU with AVX2 (and optionally AVX512) support for accelerated code paths

## Build (examples)

From the project root:

- GCC / Clang:
```bash
# AVX2 build
gcc -march=native -mavx2 -Iinclude src/*.c tests/test_add.c -o test_add
```
```bash
# AVX-512 (if supported)
gcc -march=native -mavx512f -mavx512bw -Iinclude src/*.c tests/test_add.c -o test_add_avx512
```

## Project Structure
/include      # headers
/src          # core library implementation
/tests        # small unit tests

## Contributing
Contributions are welcome but not expected, this is more of a personal project, and it kinda sucks.
Pull requests or issues are fine, but please don’t expect full support.

## License
This project is licensed under the MIT License.
Provided as-is for research/educational use.
