#ifndef _NUMMETH_H_
#define _NUMMETH_H_

#include <stdint.h>
#include <time.h>
#include <math.h>

/* ===========================================================
 * Section 0: Factorial and Complex Number Arithmetic
 * ===========================================================
 */

/**
 * @brief Compute factorial of a non-negative integer.
 *
 * @param n Input number.
 * @return n! (factorial of n)
 */
static inline unsigned long factorial(unsigned long n) {
    if (n == 0 || n == 1)
        return 1;

    unsigned long result = n;
    while (--n > 1) {
        result *= n;
    }
    return result;
}

/**
 * @brief Struct to represent complex numbers with double precision.
 */
typedef struct {
    double real;
    double imag;
} Complex;

/**
 * @brief Add two complex numbers.
 */
static inline Complex Complex_Add(Complex x, Complex y) {
    Complex result = {
        .real = x.real + y.real,
        .imag = x.imag + y.imag
    };
    return result;
}

/**
 * @brief Subtract two complex numbers.
 */
static inline Complex Complex_Sub(Complex x, Complex y) {
    Complex result = {
        .real = x.real - y.real,
        .imag = x.imag - y.imag
    };
    return result;
}

/**
 * @brief Multiply two complex numbers.
 */
static inline Complex Complex_Mul(Complex x, Complex y) {
    Complex result = {
        .real = x.real * y.real - x.imag * y.imag,
        .imag = x.real * y.imag + x.imag * y.real
    };
    return result;
}

/**
 * @brief Compute complex conjugate.
 */
static inline Complex Complex_Conj(Complex x) {
    Complex result = {
        .real = x.real,
        .imag = -x.imag
    };
    return result;
}

/**
 * @brief Compute modulus (magnitude) of a complex number.
 */
static inline double Complex_Mod(Complex x) {
    return sqrt(x.real * x.real + x.imag * x.imag);
}

/* ===========================================================
 * Section 1: Linear Congruential Generator (LCG)
 * ===========================================================
 *
 * Algorithm:
 *   X_{n+1} = (LCG_A * X_n + LCG_C) mod 2^64
 *
 * Usage:
 *   double r;
 *   uint64_t s = lcg_random(&r, 1);
 *
 *   double arr[10];
 *   s = lcg_random(arr, 10);
 *
 *   lcg_seed(12345ULL);     // manually seed
 */

#define LCG_A       6364136223846793005ULL
#define LCG_C       1442695040888963407ULL
#define LCG_M       (1ULL << 64)
#define PRE_SEED    0x9E3779B97F4A7C15ULL // 64-bit golden ratio constant

static uint64_t lcg_state = 0ULL;
static int user_seeded = 0;

/**
 * @brief Seed the LCG with a user-defined or PRE_SEED.
 *
 * @param seed User-defined seed. If 0, uses PRE_SEED.
 */
static inline void lcg_seed(uint64_t seed) {
    lcg_state = (seed != 0ULL) ? seed : PRE_SEED;
    user_seeded = 1;
}

/**
 * @brief Fill an array with random doubles in [0, 1).
 *
 * @param out Output array to fill.
 * @param n   Number of random numbers to generate.
 * @return    The seed used for generation.
 */
static inline uint64_t lcg_random(double *out, size_t n) {
    if (!user_seeded && lcg_state == 0ULL) {
        lcg_state = PRE_SEED; // Auto-seed
    }
    uint64_t used_seed = lcg_state;
    for (size_t i = 0; i < n; ++i) {
        lcg_state = (LCG_A * lcg_state + LCG_C);
        out[i] = (lcg_state >> 11) * (1.0 / (1ULL << 53));
    }
    return used_seed;
}

/**
 * @brief Test any LCG params: X_{n+1} = (A * X_n + C) mod M
 *
 * @param out Output array to fill.
 * @param n   Number of random numbers to generate.
 * @param A   Multiplier.
 * @param C   Increment.
 * @param M   Modulus.
 */
static inline uint64_t lcg_test_random(
    double *out,
    uint64_t A,
    uint64_t C,
    uint64_t M,
    uint64_t seed,
    size_t n
) {
    uint64_t state = (seed != 0ULL) ? seed : (uint64_t)PRE_SEED;
    uint64_t used_seed = state;

    for (size_t i = 0; i < n; ++i) {
        // Use bitwise AND if M is a power of two, else use modulo
        if ((M & (M - 1)) == 0 && M != 0) {
            state = (A * state + C) & (M - 1);
        } else {
            state = (A * state + C) % M;
        }
        // scale to [0, 1) even if M != 2^64
        out[i] = (double)state / (double)M;
    }
    return used_seed;
}

#endif  // _NUMMETH_H_
