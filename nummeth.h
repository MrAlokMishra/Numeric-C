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

#define LU_SUCCESS       0
#define LU_SINGULAR      -2
#define LU_INVALID       -1
#define LU_PIVOT_THRESHOLD 1e-12

/**
 * @brief LU decomposition with partial pivoting
 *
 * Decomposes square matrix A into PA = LU, using Gauss-Jordan Elemination.
 * L is unit lower-triangular, U is upper-triangular,
 * P is represented as a permutation array perm[].
 *
 * @param A Pointer to square matrix to decompose
 * @param L Pointer to pre-allocated L matrix
 * @param U Pointer to pre-allocated U matrix
 * @param perm Permutation array of size rows
 * @param num_swaps Output: number of row swaps
 * @return LU_SUCCESS on success,
 *         LU_SINGULAR if a near-zero pivot was detected,
 *         LU_INVALID if input is invalid
 */

int mat_lu_decompose(const Matrix *A, Matrix *L, Matrix *U,
                     size_t *perm, int *num_swaps) {
    if (!A || !A->data || !L || !L->data || !U || !U->data || !perm)
        return LU_INVALID;
    if (A->rows != A->cols || L->rows != L->cols || U->rows != U->cols)
        return LU_INVALID;

    size_t n = A->rows;
    int status = 0;

    if (num_swaps)
        *num_swaps = 0;

    // Initialize permutation and copy A into U
    for (size_t i = 0; i < n; i++) {
        perm[i] = i;
        for (size_t j = 0; j < n; j++) {
            U->data[i * n + j] = A->data[i * n + j];
            L->data[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    size_t pivot_row = 0, pivot_col = 0;

    while (pivot_row < n && pivot_col < n) {
        // Find max pivot in current column
        double max_val = 0.0;
        size_t max_row = pivot_row;
        for (size_t i = pivot_row; i < n; i++) {
            double val = fabs(U->data[i * n + pivot_col]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        if (max_val < LU_PIVOT_THRESHOLD) {
            // No good pivot, skip column
            pivot_col++;
            continue;
        }

        // If max_row != pivot_row
        if (max_row != pivot_row) {
            // Swap rows in U
            for (size_t j = 0; j < n; j++) {
                double tmp = U->data[pivot_row * n + j];
                U->data[pivot_row * n + j] = U->data[max_row * n + j];
                U->data[max_row * n + j] = tmp;
            }
            // Swap L columns 0..pivot_row-1
            for (size_t j = 0; j < pivot_row; j++) {
                double tmp = L->data[pivot_row * n + j];
                L->data[pivot_row * n + j] = L->data[max_row * n + j];
                L->data[max_row * n + j] = tmp;
            }
            // Swap permutation
            size_t tmp_p = perm[pivot_row];
            perm[pivot_row] = perm[max_row];
            perm[max_row] = tmp_p;

            if (num_swaps)  // only update if provided
                (*num_swaps)++;
        }

        // Eliminate below pivot
        for (size_t i = pivot_row + 1; i < n; i++) {
            double mult = U->data[i * n + pivot_col] / U->data[pivot_row * n + pivot_col];
            L->data[i * n + pivot_row] = mult;
            U->data[i * n + pivot_col] = 0.0; // explicitly zero lower part
            for (size_t j = pivot_col + 1; j < n; j++) {
                U->data[i * n + j] -= mult * U->data[pivot_row * n + j];
            }
        }

        pivot_row++;
        pivot_col++;
    }

    // Check if any diagonal of U is too small
    for (size_t i = 0; i < n; i++) {
        if (fabs(U->data[i * n + i]) < LU_PIVOT_THRESHOLD)
            status = LU_SINGULAR;
    }

    return status;
}


/**
 * @brief Apply a row permutation to a matrix.
 *
 * Copies rows from `src` to `dest` according to `perm`.
 *
 * @param src Source matrix (rows x cols).
 * @param dest Destination matrix (pre-allocated, same size as src).
 * @param perm Permutation array of length src->rows.
 * @return 0 on success, -1 if pointers or dimensions are invalid.
 */
static inline int mat_apply_perm(const Matrix *src, Matrix *dest, const size_t *perm) {
    if (!src || !src->data || !dest || !dest->data || !perm) return -1;
    if (src->rows != dest->rows || src->cols != dest->cols) return -1;

    for (size_t i = 0; i < src->rows; i++)
        for (size_t j = 0; j < src->cols; j++)
            dest->data[i * src->cols + j] = src->data[perm[i] * src->cols + j];

    return 0;
}



/**
 * @brief Solve a linear system Ax = b using LU decomposition (with permutation).
 *
 * This function performs the LU decomposition of A internally.
 *
 * @param A Coefficient matrix (n x n)
 * @param b Right-hand side vector (length n)
 * @param x Solution vector (pre-allocated, length n)
 * @return LU_SUCCESS on success,
 *         LU_SINGULAR if U has a near-zero pivot,
 *         LU_INVALID if inputs are invalid
 */
int mat_lu_solve(const Matrix *A, const double *b, double *x) {
    if (!A || !b || !x) return LU_INVALID;
    if (A->rows != A->cols) return LU_INVALID;

    size_t n = A->rows;

    // Allocate L, U, and permutation
    Matrix *L = mat_new(n, n);
    Matrix *U = mat_new(n, n);
    size_t *perm = malloc(n * sizeof *perm);

    if (!L || !U || !perm) {
        if (L) mat_free(L);
        if (U) mat_free(U);
        if (perm) free(perm);
        return LU_INVALID;
    }

    // Perform LU decomposition
    int status = mat_lu_decompose(A, L, U, perm , NULL);
    if (status != LU_SUCCESS) {
        mat_free(L);
        mat_free(U);
        free(perm);
        return status; // could be LU_SINGULAR or LU_INVALID
    }

    // Step 1: Apply permutation to b → b'
    double *bp = malloc(n * sizeof *bp);
    if (!bp) {
        mat_free(L);
        mat_free(U);
        free(perm);
        return LU_INVALID;
    }
    for (size_t i = 0; i < n; i++)
        bp[i] = b[perm[i]];

    // Step 2: Forward substitution Ly = b'
    double *y = malloc(n * sizeof *y);
    if (!y) {
        free(bp);
        mat_free(L);
        mat_free(U);
        free(perm);
        return LU_INVALID;
    }

    for (size_t i = 0; i < n; i++) {
        double sum = bp[i];
        for (size_t j = 0; j < i; j++) {
            sum -= L->data[i * n + j] * y[j];
        }
        y[i] = sum;  // L has unit diagonal
    }

    // Step 3: Backward substitution Ux = y
    for (size_t i = n - 1; i < n; i--) {
        double sum = y[i];
        for (size_t j = i + 1; j < n; j++) {
            sum -= U->data[i * n + j] * x[j];
        }
        if (fabs(U->data[i * n + i]) < LU_PIVOT_THRESHOLD) {
            free(bp); free(y);
            mat_free(L); mat_free(U); free(perm);
            return LU_SINGULAR;
        }
        x[i] = sum / U->data[i * n + i];
    }

    // Cleanup
    free(bp);
    free(y);
    mat_free(L);
    mat_free(U);
    free(perm);

    return LU_SUCCESS;
}


/**
 * @brief Compute determinant of a square matrix using LU decomposition
 *
 * This function performs LU decomposition of A and computes the determinant as
 * det(A) = (-1)^num_swaps * product of diagonal elements of U.
 *
 * @param A Pointer to square matrix (Matrix*) created with mat_new()
 * @return Determinant value as double:
 *         - Returns 0.0 if the matrix is singular.
 *         - Returns NAN if A is NULL, non-square, or memory allocation fails.
 */
double mat_determinant(const Matrix *A) {
    if (!A || !A->data || A->rows != A->cols)
        return NAN;

    size_t n = A->rows;
    
    Matrix *L = mat_new(n, n);
    Matrix *U = mat_new(n, n);
    size_t *perm = malloc(n * sizeof(size_t));
    int num_swaps = 0;

    if (!L || !U || !perm) {
        if (L) mat_free(L);
        if (U) mat_free(U);
        free(perm);
        return NAN;  // memory allocation error
    }

    int status = mat_lu_decompose(A, L, U, perm, &num_swaps);
    if (status == LU_SINGULAR) {
        mat_free(L);
        mat_free(U);
        free(perm);
        return 0.0;  // singular matrix
    } else if (status != LU_SUCCESS) {
        mat_free(L);
        mat_free(U);
        free(perm);
        return NAN;  // invalid input or other error
    }

    // Compute det(U)
    double detU = 1.0;
    for (size_t i = 0; i < n; i++)
        detU *= U->data[i * n + i];

    // determinant = (-1)^num_swaps * detU
    double det = (num_swaps % 2 == 0) ? detU : -detU;

    // Cleanup
    mat_free(L);
    mat_free(U);
    free(perm);

    return det;
}




/**
 * @brief Compute the inverse of a square matrix A and store it in inv.
 *
 * The function computes the inverse by solving n linear systems A * x = e_i for
 * each canonical basis vector e_i (i = 0..n-1) and placing each solution vector
 * x as the i-th column of inv. Temporary working arrays are allocated internally.
 *
 * @param  A   Pointer to the input Matrix to invert.
 * @param  inv Pointer to the Matrix that will receive the inverse.
 *
 * @return LU_SUCCESS on success.
 *         LU_INVALID If any pointer is NULL, A is not square, inv has mismatched
 *                    dimensions, or an internal memory allocation for temporary
 *                    buffers fails.
 */
int mat_inverse(const Matrix *A, Matrix *inv) {
    if (!A || !A->data || !inv || !inv->data) return LU_INVALID;
    if (A->rows != A->cols) return LU_INVALID;
    if (inv->rows != A->rows || inv->cols != A->cols) return LU_INVALID;

    size_t n = A->rows;
    double *b = malloc(n * sizeof *b);
    double *x = malloc(n * sizeof *x);
    if (!b || !x) {
        if(b) free(b);
        if(x) free(x);
        return LU_INVALID;
    }

    for (size_t col = 0; col < n; ++col) {
        for (size_t i = 0; i < n; ++i) b[i] = (i == col) ? 1.0 : 0.0;

        int status = mat_lu_solve(A, b, x);
        if (status != LU_SUCCESS) {
            if(b) free(b);
            if(x) free(x);
            return status;
        }

        for (size_t i = 0; i < n; ++i)
            inv->data[i * n + col] = x[i];
    }

    free(b);
    free(x);
    return LU_SUCCESS;
}



#endif  // _NUMMETH_H_
