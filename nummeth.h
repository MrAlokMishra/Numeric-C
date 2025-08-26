#ifndef _NUMMETH_H_
#define _NUMMETH_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/**
 * @section Section 0: Factorial and Complex Number Arithmetic
 *
 * @brief Basic arithmetic utilities for integers and complex numbers.
 *
 * @details
 * Provides:
 * - Factorial computation for non-negative integers.
 * - Complex number operations: addition, subtraction, multiplication, conjugate, and modulus.
 */



/**
 * @brief Compute the factorial of a non-negative integer.
 *
 * @param n Non-negative integer input.
 * @return Factorial of n (n!) as an unsigned long, or (unsigned long)-1 on invalid input.
 *
 * @note Returns -1 if n is negative; valid input should be >= 0.
 */
static inline unsigned long factorial(unsigned long n) {
    
    if (n < 0)
        return -1 ;
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

/**
 * @section Section 1:  LCG Linear Congruential Generator (LCG)
 * @brief Implements a 64-bit Linear Congruential Generator for pseudo-random number generation.
 *
 * @details
 * The LCG algorithm generates a sequence of pseudo-random numbers using the recurrence relation:
 *   X_{n+1} = (LCG_A * X_n + LCG_C) mod 2^64
 *
 * - @c lcg_random() fills a double or an array of doubles with random values in [0, 1).
 * - @c lcg_seed() sets the initial seed for the generator.
 *
 * @usage
 * @code
 * double r;
 * uint64_t s = lcg_random(&r, 1); // Generate one random double
 *
 * double arr[10];
 * s = lcg_random(arr, 10);        // Generate an array of random doubles
 *
 * lcg_seed(12345ULL);             // Manually set the seed
 * @endcode
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
 * @param out Output array to fill (must not be NULL).
 * @param n   Number of random numbers to generate.
 * @return    The seed used for generation.
 *            Returns 0ULL if @p out is NULL.
 */
static inline uint64_t lcg_random(double *out, size_t n) {
    if (!out) return 0ULL;

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
 * @param out  Output array to fill (must not be NULL).
 * @param A    Multiplier.
 * @param C    Increment.
 * @param M    Modulus.
 * @param seed Initial seed (0ULL uses default PRE_SEED).
 * @param n    Number of random numbers to generate.
 * @return     The seed used for generation.
 *             Returns 0ULL if @p out is NULL.
 */
static inline uint64_t lcg_test_random(
    double *out,
    uint64_t A,
    uint64_t C,
    uint64_t M,
    uint64_t seed,
    size_t n
) {
    if (!out) return 0ULL;

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


/**
 * @brief Generate samples from an arbitrary distribution using inverse CDF sampling.
 *
 * @param inv_cdf  Function pointer to inverse CDF: double f(double u).
 * @param out      Output array of samples.
 * @param n        Number of samples to generate.
 * @return         The seed used for generation on success,
 *                 or 0ULL if an error occurs (e.g., NULL pointer or malloc failure).
 *
 */
static inline uint64_t lcg_inverse_sample(double (*inv_cdf)(double), double *out, size_t n) {
    if (!inv_cdf || !out) return 0ULL;

    // Step 1: generate uniform(0,1) numbers
    double *uniforms = malloc(n * sizeof(double));
    if (!uniforms) return 0ULL;

    uint64_t seed_used = lcg_random(uniforms, n);

    // Step 2: apply inverse CDF to each uniform
    for (size_t i = 0; i < n; ++i) {
        out[i] = inv_cdf(uniforms[i]);
    }

    free(uniforms);
    return seed_used;
}


/**
 * @section Section 2: Matrix Struct and Basic Operations
 *
 * Defines a dense matrix type with row-major storage and basic operations.
 *
 * @details
 * Includes memory management (`mat_new`, `mat_free`), I/O (`mat_load`, `mat_save`, `mat_read`, `mat_print`),
 * element-wise arithmetic (`mat_add`, `mat_sub`), and matrix multiplication (`mat_mul`).
 * All functions check pointers and dimensions, returning error codes for invalid input.
 * Provides the foundation for LU decomposition, solving systems, determinants, and inversion.
 */



/**
 * @brief Dense matrix in row-major order.
 */
typedef struct {
    size_t rows;   /**< number of rows */
    size_t cols;   /**< number of columns */
    double *data;  /**< pointer to row-major data (rows * cols) */
} Matrix;

/**
 * @brief Create a new matrix (allocated).
 *
 * @param rows Number of rows (must be > 0).
 * @param cols Number of columns (must be > 0).
 * @return Pointer to allocated Matrix or NULL on failure.
 */
static inline Matrix *mat_new(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0)
        return NULL;

    Matrix *m = (Matrix *)malloc(sizeof(Matrix));
    if (!m)
        return NULL;

    m->rows = rows;
    m->cols = cols;
    m->data = (double *)calloc(rows * cols, sizeof(double));
    if (!m->data) {
        free(m);
        return NULL;
    }
    return m;
}

/**
 * @brief Free a matrix previously created with mat_new.
 *
 * @param m Pointer to Matrix to free (may be NULL).
 */
static inline void mat_free(Matrix *m) {
    if (m) {
        free(m->data);
        m->data = NULL;
        free(m);
    }
}

/**
 * @brief Load matrix values from a file into a pre-created matrix.
 *
 * @param filename Path to input file.
 * @param m       Pointer to pre-allocated Matrix (m->data must be allocated).
 * @return 0 on success,
 *         -1 if filename, m, or m->data is NULL,
 *         -2 if file cannot be opened,
 *         -3 if not enough numbers in file,
 *         -4 if extra numbers in file.
 */
static inline int mat_load(const char *filename, Matrix *m) {
    if (!filename || !m || !m->data) return -1;

    FILE *file = fopen(filename, "r");
    if (!file) {
        return -2;
    }

    for (size_t i = 0; i < m->rows; i++) {
        for (size_t j = 0; j < m->cols; j++) {
            if (fscanf(file, "%lf", &m->data[i * m->cols + j]) != 1) {
                fclose(file);
                return -3;
            }
        }
    }

    double tmp;
    if (fscanf(file, "%lf", &tmp) == 1) {
        fclose(file);
        return -4;
    }

    fclose(file);
    return 0;
}

/**
 * @brief Read matrix values from standard input into a pre-created matrix.
 *
 * @param m Pointer to pre-allocated Matrix (m->data must be allocated).
 * @return 0 on success,
 *        -1 if m or m->data is NULL,
 *        -2 if invalid input.
 */
static inline int mat_read(Matrix *m) {
    if (!m || !m->data) return -1;

    for (size_t i = 0; i < m->rows; i++) {
        for (size_t j = 0; j < m->cols; j++) {
            if (scanf("%lf", &m->data[i * m->cols + j]) != 1) {
                return -2;
            }
        }
    }
    return 0;
}

/**
 * @brief Save matrix to a file in row-major order.
 *
 * @param filename Path to output file.
 * @param m       Pointer to Matrix to save.
 * @return 0 on success,
 *        -1 if filename, m, or m->data is NULL,
 *        -2 if file cannot be opened for writing.
 */
static inline int mat_save(const char *filename, const Matrix *m) {
    if (!filename || !m || !m->data) return -1;

    FILE *file = fopen(filename, "w");
    if (!file) {
        return -2;
    }

    for (size_t i = 0; i < m->rows; i++) {
        for (size_t j = 0; j < m->cols; j++) {
            fprintf(file, "%lf ", m->data[i * m->cols + j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
    return 0;
}

/**
 * @brief Print matrix to standard output.
 *
 * @param m Pointer to Matrix to print.
 */
static inline void mat_print(const Matrix *m) {
    if (!m || !m->data) return;

    for (size_t i = 0; i < m->rows; i++) {
        for (size_t j = 0; j < m->cols; j++) {
            printf("%lf ", m->data[i * m->cols + j]);
        }
        printf("\n");
    }
}

/**
 * @brief Add two matrices element-wise.
 * 
 * @param a Pointer to the first matrix.
 * @param b Pointer to the second matrix.
 * @param result Pointer to the pre-created matrix to store the result.
 * @return 0 on success,
 *        -1 if dimensions do not match,
 *        -2 if any input pointer (a, b, result or their data pointers) is NULL.
 */
int mat_add(const Matrix *a, const Matrix *b, Matrix *result) {
    if (!a || !a->data || !b || !b->data || !result || !result->data) return -2;
    if (a->rows != b->rows || a->cols != b->cols) return -1;
    if (result->rows != a->rows || result->cols != a->cols) return -1;

    for (size_t i = 0; i < a->rows; i++) {
        for (size_t j = 0; j < a->cols; j++) {
            result->data[i * a->cols + j] = a->data[i * a->cols + j] + b->data[i * b->cols + j];
        }
    }
    return 0;
}

/**
 * @brief Subtract two matrices element-wise (a - b).
 * 
 * @param a Pointer to the first matrix.
 * @param b Pointer to the second matrix.
 * @param result Pointer to the pre-created matrix to store the result.
 * @return 0 on success,
 *        -1 if dimensions do not match,
 *        -2 if any input pointer (a, b, result or their data pointers) is NULL.
 */
int mat_sub(const Matrix *a, const Matrix *b, Matrix *result) {
    if (!a || !a->data || !b || !b->data || !result || !result->data) return -2;
    if (a->rows != b->rows || a->cols != b->cols) return -1;
    if (result->rows != a->rows || result->cols != a->cols) return -1;

    for (size_t i = 0; i < a->rows; i++) {
        for (size_t j = 0; j < a->cols; j++) {
            result->data[i * a->cols + j] = a->data[i * a->cols + j] - b->data[i * b->cols + j];
        }
    }
    return 0;
}

/**
 * @brief Multiply two matrices (matrix product).
 * 
 * @param a Pointer to the first matrix.
 * @param b Pointer to the second matrix.
 * @param result Pointer to the pre-created matrix to store the result.
 * @return 0 on success,
 *        -1 if inner dimensions do not match or result has wrong size,
 *        -2 if any input pointer (a, b, result or their data pointers) is NULL.
 */
int mat_mul(const Matrix *a, const Matrix *b, Matrix *result) {
    if (!a || !a->data || !b || !b->data || !result || !result->data) return -2;
    if (a->cols != b->rows) return -1;
    if (result->rows != a->rows || result->cols != b->cols) return -1;

    for (size_t i = 0; i < a->rows; i++) {
        for (size_t j = 0; j < b->cols; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < a->cols; k++) {
                sum += a->data[i * a->cols + k] * b->data[k * b->cols + j];
            }
            result->data[i * b->cols + j] = sum;
        }
    }
    return 0;
}




/**
 * @section Section 3: Matrix Application functions
/**
 * @section Section 3: Matrix Application functions
 * @brief LU, permutation, solve, determinant and inverse utilities for Matrix.
 *
 * @details
 * - mat_lu_decompose: PA = LU with partial pivoting (perm[], num_swaps).
 * - mat_apply_perm: apply row permutation to a matrix.
 * - mat_lu_solve: solve Ax = b using LU (allocates temporaries).
 * - mat_determinant: det(A) = (-1)^swaps * prod(diag(U)); 0 for singular, NAN on error.
 * - mat_inverse: compute A^{-1} by solving for canonical basis vectors.
 *
 * Inputs are validated; functions return LU_SUCCESS, LU_SINGULAR, or LU_INVALID.
 */

#define LU_SUCCESS        0
#define LU_INVALID       -1
#define LU_SINGULAR      -2
#define LU_PIVOT_THRESHOLD 1e-12

/**
 * @brief In-place LU decomposition with partial pivoting.
 *
 * Factorizes a square matrix A into PA = LU, where:
 * - L is unit lower-triangular (1s on diagonal, multipliers stored below diagonal in A),
 * - U is upper-triangular (stored in and above the diagonal in A),
 * - P is a row permutation matrix represented by `perm[]`.
 *
 * The decomposition is performed in-place: the input matrix A is overwritten
 * with the combined LU factors. No additional L or U matrices are allocated.
 *
 * @param A          Pointer to a square matrix (modified in-place to contain LU).
 * @param perm       Output array of length A->rows, representing row permutations.
 * @param num_swaps  Output: number of row swaps performed (may be NULL if unused).
 *
 * @return LU_SUCCESS  on success,
 *         LU_SINGULAR if a near-zero pivot was detected,
 *         LU_INVALID  if input is invalid (null pointers or non-square matrix).
 */
int mat_lu_decompose(Matrix *A, size_t *perm, int *num_swaps) {
    if (!A || !A->data || !perm) return LU_INVALID;
    if (A->rows != A->cols) return LU_INVALID;

    size_t n = A->rows;
    if (num_swaps) *num_swaps = 0;

    // initialize permutation
    for (size_t i = 0; i < n; i++)
        perm[i] = i;

    for (size_t k = 0; k < n; k++) {
        // --- Pivot selection ---
        size_t max_row = k;
        double max_val = fabs(A->data[k * n + k]);
        for (size_t i = k + 1; i < n; i++) {
            double val = fabs(A->data[i * n + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        if (max_val < LU_PIVOT_THRESHOLD)
            return LU_SINGULAR;

        // --- Row swap if necessary ---
        if (max_row != k) {
            for (size_t j = 0; j < n; j++) {
                double tmp = A->data[k * n + j];
                A->data[k * n + j] = A->data[max_row * n + j];
                A->data[max_row * n + j] = tmp;
            }
            size_t tmp_p = perm[k];
            perm[k] = perm[max_row];
            perm[max_row] = tmp_p;

            if (num_swaps) (*num_swaps)++;
        }

        // --- Elimination step (store multipliers in lower triangle) ---
        for (size_t i = k + 1; i < n; i++) {
            double mult = A->data[i * n + k] / A->data[k * n + k];
            A->data[i * n + k] = mult;
            for (size_t j = k + 1; j < n; j++) {
                A->data[i * n + j] -= mult * A->data[k * n + j];
            }
        }
    }

    return LU_SUCCESS;
}




/**
 * @brief Solve Ax = b in-place using LU decomposition with partial pivoting.
 *
 * On return:
 *   - A is overwritten by its LU decomposition
 *   - b is overwritten with the solution x
 *
 * @param A (n x n) coefficient matrix, modified to contain LU
 * @param b (length n) RHS vector, overwritten with solution x
 * @return LU_SUCCESS, LU_SINGULAR, or LU_INVALID
 */
int mat_lu_solve(Matrix *A, double *b) {
    if (!A || !b) return LU_INVALID;
    if (A->rows != A->cols) return LU_INVALID;

    size_t n = A->rows;
    size_t *perm = malloc(n * sizeof *perm);
    if (!perm) return LU_INVALID;

    // Perform in-place LU decomposition
    int status = mat_lu_decompose(A, perm, NULL);
    if (status != LU_SUCCESS) {
        free(perm);
        return status;
    }

    // Apply row permutations to b directly
    for (size_t i = 0; i < n; i++) {
        if (perm[i] != i) {
            double tmp = b[i];
            b[i] = b[perm[i]];
            b[perm[i]] = tmp;
        }
    }

   // Forward substitution: solve L y = P b
for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < i; j++) {   // strictly lower part of A
        b[i] -= A->data[i * n + j] * b[j];
    }
    // L has implicit 1s on diag → nothing else needed
}

// Backward substitution: solve U x = y
for (long i = (long)n - 1; i >= 0; i--) {
    for (size_t j = i + 1; j < n; j++) {   // strictly upper part of A
        b[i] -= A->data[i * n + j] * b[j];
    }
    double diag = A->data[i * n + i];   // U's diagonal
    if (fabs(diag) < LU_PIVOT_THRESHOLD) {
        free(perm);
        return LU_SINGULAR;
    }
    b[i] /= diag;
}
    free(perm);
    return LU_SUCCESS;
}


/**
 * @brief Compute determinant of a square matrix using in-place LU decomposition
 *
 * This function performs LU decomposition of A (in-place) and computes:
 *      det(A) = (-1)^num_swaps * product of diagonal elements of U.
 *
 * @param A Pointer to square matrix (Matrix*) created with mat_new()
 * @return Determinant value as double:
 *         - Returns 0.0 if the matrix is singular.
 *         - Returns NAN if A is NULL, non-square, or memory allocation fails.
 *
 * Note: This function modifies A.
 */
static double mat_determinant(Matrix *A) {
    if (!A || !A->data || A->rows != A->cols)
        return NAN;
    if (A->rows != A->cols) return NAN; // non-square

    size_t n = A->rows;
    size_t *perm = malloc(n * sizeof(size_t));
    if (!perm) return NAN;

    int num_swaps = 0;
    int status = mat_lu_decompose(A, perm, &num_swaps);
    if (status == LU_SINGULAR) {
        free(perm);
        return 0.0; // Singular matrix → det = 0
    } else if (status != LU_SUCCESS) {
        free(perm);
        return NAN; // Error
    }

    // determinant = (-1)^num_swaps * product of diagonal of U
    double det = (num_swaps % 2 == 0) ? 1.0 : -1.0;
    for (size_t i = 0; i < n; i++)
        det *= A->data[i * n + i];

    free(perm);
    return det;
}



/**
 * @brief Computes the inverse of a square matrix using LU decomposition.
 *
 * This function performs an in-place LU decomposition of the input matrix A,
 * then explicitly computes the inverse of the upper-triangular matrix U,
 * and solves for A⁻¹ by back-substitution and applying the permutation.
 * 
 * The algorithm follows the LAPACK-style approach:
 *   A⁻¹ P⁻¹ L = U⁻¹
 *
 * @param[in,out] A   Pointer to the input matrix (will be LU-decomposed in place).
 * @param[out]    inv Pointer to the output matrix where the inverse is stored.
 *
 * @return LU_SUCCESS if inversion succeeds,
 *         LU_INVALID if input is invalid or allocation fails,
 *         other LU_* error codes from mat_lu_decompose.
 *
 * Note : Function changes A
 */
int mat_inverse(Matrix *A, Matrix *inv) {
    if (!A || !A->data || !inv || !inv->data) return LU_INVALID;
    if (A->rows != A->cols) return LU_INVALID;

    size_t n = A->rows;
    if (inv->rows != n || inv->cols != n) return LU_INVALID;

    // Allocate permutation array
    size_t *perm = malloc(n * sizeof *perm);
    if (!perm) return LU_INVALID;

    // Step 1: LU decomposition in-place
    int status = mat_lu_decompose(A, perm, NULL);
    if (status != LU_SUCCESS) {
        free(perm);
        return status;
    }


    // Step 2: Initialize Uinv as identity
    Matrix *Uinv = mat_new(n, n);
    if (!Uinv->data) {
        free(perm);
        return LU_INVALID;
    }

    for (size_t i = 0; i < n * n; i++)
        Uinv->data[i] = 0.0;

    for (size_t i = 0; i < n; i++)
        Uinv->data[i * n + i] = 1.0;

    // Step 3: Compute U^-1 (upper-triangular inversion)
    for (long i = (long)n - 1; i >= 0; i--) {
        Uinv->data[i * n + i] = 1 / A->data[i * n + i]; // diagonal

        for (long j = i - 1; j >= 0; j--) {
            double sum = 0.0;

            for (long k = j + 1; k <= i; k++)
                sum += A->data[j * n + k] * Uinv->data[k * n + i];

            Uinv->data[j * n + i] = -sum / A->data[j * n + j];
        }
    }


    // Step 4: Forward substitution X L = U^-1
    for (size_t row = 0; row < n; row++) {             
        for (long col = (long)n - 2; col >= 0; col--) { 
            double sum = 0.0;

            for (size_t k = col + 1; k < n; k++)
                sum += A->data[k * n + col] * Uinv->data[row * n + k];

            Uinv->data[row * n + col] -= sum;
        }
    }

    // Step 5: Apply permutation: inv = X P
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            inv->data[i * n + j] = Uinv->data[i * n + perm[j]];

    mat_free(Uinv);
    free(perm);

    return LU_SUCCESS;
}


/**
 * @brief In-place Cholesky decomposition of a Hermitian positive definite matrix.
 *
 * @param A Matrix to decompose (must be square and strictly Hermitian).
 *          On success, the lower triangle contains L such that A = L L^T,
 *          and the upper triangle is zeroed out.
 *
 * @note The input must be strictly Hermitian and positive definite.
 *       If the matrix is not positive definite, decomposition fails.
 *
 * @return 0 on success, or k (1-based index) if the leading minor of order k
 *         is not positive definite.
 */

static int mat_chol_hpd(Matrix *A) {
    if (!A || A->rows != A->cols) return -1;

    size_t n = A->rows;

    for (size_t j = 0; j < n; j++) {
        double sum = A->data[j * n + j];

        // subtract L[j,0..j-1]^2
        for (size_t k = 0; k < j; k++) {
            double l_jk = A->data[j * n + k];
            sum -= l_jk * l_jk;
        }

        if (sum <= 0.0) {
            return (int)(j + 1); // not positive definite
        }

        double l_jj = sqrt(sum);
        A->data[j * n + j] = l_jj;

        // update column j below the diagonal
        for (size_t i = j + 1; i < n; i++) {
            double s = A->data[i * n + j];
            for (size_t k = 0; k < j; k++) {
                s -= A->data[i * n + k] * A->data[j * n + k];
            }
            A->data[i * n + j] = s / l_jj;
        }
    }

    // optional: zero out upper triangle
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            A->data[i * n + j] = 0.0;
        }
    }

    return 0;
}

/**
 * @brief Solve A x = b using Cholesky factorization (in place).
 *
 * @param A Matrix to decompose (must be square and strictly Hermitian 
 *          positive definite). Will be overwritten by its Cholesky factor L.
 * @param b Right-hand side vector, overwritten with the solution x.
 *
 * @note The input matrix must be strictly Hermitian and positive definite.
 *       If factorization fails, the system cannot be solved.
 *
 * @return 0 on success, or k (1-based index) if the leading minor of order k
 *         is not positive definite (factorization failed).
 */
static int mat_chol_solve(Matrix *A, double *b) {
    if (!A || A->rows != A->cols || !b || !A->data) return -1;

    size_t n = A->rows;

    // Step 1: Cholesky decomposition in-place
    int status = mat_chol_hpd(A);
    if (status != 0) {
        return status; // not positive definite
    }

    // Step 2: Forward substitution (L y = b)
    for (size_t i = 0; i < n; i++) {
        double sum = b[i];
        for (size_t j = 0; j < i; j++) {
            sum -= A->data[i * n + j] * b[j];
        }
        b[i] = sum / A->data[i * n + i];
    }

    // Step 3: Back substitution (L^T x = y)
    for (long i = (long)n - 1; i >= 0; i--) {
        double sum = b[i];
        for (size_t j = i + 1; j < n; j++) {
            sum -= A->data[j * n + i] * b[j];
        }
        b[i] = sum / A->data[i * n + i];
    }

    return 0;
}

/**
 * @brief Compute determinant of a Hermitian positive definite matrix using Cholesky factorization.
 *
 * @param A Matrix (must be square and strictly Hermitian positive definite).
 *          On success, A is overwritten by its Cholesky factor L.
 *
 * @return Determinant of A on success, or NAN if factorization fails.
 */
double mat_det_hpd(Matrix* A) {
    if (!A || A->rows != A->cols || !A->data) return -1;
    size_t n = A->rows;
    int status = mat_chol_hpd(A);
    if (status != 0) {
        return NAN;  // failed factorization -> not HPD
    }

    double prod = 1.0;
    for (size_t i = 0; i < n; i++) {
        double Lii = A->data[i*n + i];
        prod *= Lii;
    }
    return prod * prod;  // det(A) = (∏ diag(L))^2
}




#endif  // _NUMMETH_H_
