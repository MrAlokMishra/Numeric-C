#ifndef _NUMMETH_H_
#define _NUMMETH_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// ==================== Scalar Abstraction Layer ====================
// Switch between double and Complex by changing one typedef/define


/**
 * @section Section 0: Complex Number Arithmetic and Polynomials
 *
 * @brief Basic arithmetic utilities for integers and complex numbers.
 *
 * @details
 * Provides:
 * - Complex number operations: addition, subtraction, multiplication, conjugate, and modulus.
 */

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
 * @brief Divide two complex numbers.
 */
static inline Complex Complex_Div(Complex x, Complex y) {
    double denom = y.real*y.real + y.imag*y.imag;
    Complex result = {
        .real = (x.real*y.real + x.imag*y.imag) / denom,
        .imag = (x.imag*y.real - x.real*y.imag) / denom
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
 * @brief Compute the square root of a complex number.
 *
 * @details
 * Uses the standard formula:
 *   sqrt(x + iy) = ±(u + iv),
 * where
 *   u = sqrt((|z| + x)/2),
 *   v = sign(y) * sqrt((|z| - x)/2).
 *
 * @param z Input complex number.
 * @return Principal square root of z.
 */
static inline Complex Complex_Sqrt(Complex z) {
    Complex result;
    double r = sqrt(z.real * z.real + z.imag * z.imag); // |z|

    double u = sqrt((r + z.real) / 2.0);
    double v = sqrt((r - z.real) / 2.0);

    if (z.imag < 0)
        v = -v;

    result.real = u;
    result.imag = v;
    return result;
}

#ifdef USE_COMPLEX
    typedef Complex SCALAR;

    #define SCALAR_ZERO    (Complex){0.0, 0.0}
    #define SCALAR_ONE     (Complex){1.0, 0.0}

    #define SCALAR_ADD(a,b) Complex_Add(a,b)
    #define SCALAR_SUB(a,b) Complex_Sub(a,b)
    #define SCALAR_MUL(a,b) Complex_Mul(a,b)
    #define SCALAR_CONJ(a)  Complex_Conj(a)
    #define SCALAR_ABS(a)   Complex_Mod(a)
    #define SCALAR_DIV(a,b) Complex_Div(a,b)
    #define SCALAR_SQRT(a)  Complex_Sqrt(a)
    #define CAST_SCALAR(x) (Complex){(double)(x), 0.0}
    #define SCALAR_PRINT(x) printf("(%f %+fi)", (x).real, (x).imag)

#else
    typedef double SCALAR;

    #define SCALAR_ZERO    0.0
    #define SCALAR_ONE     1.0

    #define SCALAR_ADD(a,b) ((a)+(b))
    #define SCALAR_SUB(a,b) ((a)-(b))
    #define SCALAR_MUL(a,b) ((a)*(b))
    #define SCALAR_DIV(a,b) ((a)/(b))
    #define SCALAR_CONJ(a)  (a)
    #define SCALAR_ABS(a)   fabs(a)
    #define SCALAR_SQRT(a)  sqrt(a)
    #define CAST_SCALAR(x) ((double)(x))
    #define SCALAR_PRINT(x) printf("%f", (x))
#endif


// =====================================================
// Polynomial struct (generic SCALAR: double or Complex)
// =====================================================

/**
 * @brief Polynomial object with generic SCALAR coefficients.
 *
 * @details
 * Represents p(x) = Σ a_i x^i with degree = deg.
 * Coefficients are stored in coef[i] for x^i.
 */
typedef struct {
    int deg;        ///< Degree of polynomial
    SCALAR *coef;   ///< Coefficients array (length = deg+1)
} Polynomial;

// ---------- Init / Free ----------

/**
 * @brief Allocate and initialize a polynomial of degree n.
 * @param n Degree of polynomial.
 * @return Pointer to allocated polynomial (all coeffs = 0).
 */
Polynomial* poly_init(int n) {
    Polynomial *p = malloc(sizeof(Polynomial));
    if (!p) return NULL;
    p->deg = n;
    p->coef = calloc(n + 1, sizeof(SCALAR));
    return p;
}

/**
 * @brief Free a polynomial object.
 * @param p Polynomial to free.
 */
void poly_free(Polynomial *p) {
    if (p) {
        free(p->coef);
        free(p);
    }
}

// ---------- Set / Get ----------

/**
 * @brief Set coefficient of x^i.
 * @param p Polynomial.
 * @param i Exponent index.
 * @param value New coefficient value.
 */
void poly_set(Polynomial *p, int i, SCALAR value) {
    if (p && i >= 0 && i <= p->deg) p->coef[i] = value;
}

/**
 * @brief Get coefficient of x^i.
 * @param p Polynomial.
 * @param i Exponent index.
 * @return Coefficient value, or SCALAR_ZERO if out of range.
 */
SCALAR poly_get(const Polynomial *p, int i) {
    if (!p || i < 0 || i > p->deg) return SCALAR_ZERO;
    return p->coef[i];
}

// ---------- Evaluation ----------

/**
 * @brief Evaluate polynomial at point x using Horner's method.
 * @param p Polynomial.
 * @param x Point to evaluate.
 * @return Value p(x).
 */
SCALAR poly_eval(const Polynomial *p, SCALAR x) {
    SCALAR result = SCALAR_ZERO;
    for (int i = p->deg; i >= 0; i--) {
        result = SCALAR_ADD(SCALAR_MUL(result, x), p->coef[i]);
    }
    return result;
}

// ---------- Derivative ----------

/**
 * @brief Compute the derivative of a polynomial.
 * @param p Polynomial.
 * @return New polynomial representing p'(x).
 */
Polynomial* poly_derivative(const Polynomial *p) {
    if (p->deg <= 0) return poly_init(0);
    Polynomial *dp = poly_init(p->deg - 1);
    for (int i = 1; i <= p->deg; i++) {
        dp->coef[i - 1] = SCALAR_MUL(CAST_SCALAR(i), p->coef[i]);
    }
    return dp;
}

// ---------- Deflation ----------

/**
 * @brief Deflate polynomial by dividing by (x - root).
 * @param p Polynomial.
 * @param root Root to remove.
 * @return New polynomial of degree p->deg - 1.
 */
Polynomial* poly_deflate(const Polynomial *p, SCALAR root) {
    if (p->deg <= 0) return NULL;
    Polynomial *q = poly_init(p->deg - 1);
    SCALAR b = p->coef[p->deg];
    q->coef[p->deg - 1] = b;

    for (int i = p->deg - 1; i > 0; i--) {
        b = SCALAR_ADD(p->coef[i], SCALAR_MUL(root, b));
        q->coef[i - 1] = b;
    }
    return q;
}

// ---------- Debug Print ----------

/**
 * @brief Print polynomial to stdout.
 * @param p Polynomial.
 */
void poly_print(const Polynomial *p) {
    for (int i = p->deg; i >= 0; i--) {
        printf("(");
#ifdef USE_COMPLEX
        printf("%g + %gi", p->coef[i].real, p->coef[i].imag);
#else
        printf("%g", p->coef[i]);
#endif
        printf(")x^%d", i);
        if (i > 0) printf(" + ");
    }
    printf("\n");
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
 * @brief Generate samples from a standard normal distribution (mean=0, stddev=1).
 *        Uses Box Muller Transformation with Von-Nuemann Rejection sampling and Marsagila Trick
 * @param out Output array of samples.
 * @param n Number of samples to generate.
 *
 * @return The seed used for generation.
 */
static inline uint64_t lcg_std_gauss_sample(double *out, size_t n) {
    if (!out) return 0ULL;

    uint64_t seed_used = lcg_random(out, n);
    for (size_t i = 0; i < n; ++i) {
        // Box Muller with Von-Nuemann and Marsagila 
        double u1, u2, s;
        while (1) {
            u1 = out[i];
            lcg_random(&u2, 1);
            s = u1 * u1 + u2 * u2;
            if (s < 1.0 && s != 0.0) {
            break;
            }
        }
        out[i] = u1 * sqrt(-2.0 * log(s) / s);
    }
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
    if (!file) return -2;

    for (size_t i = 0; i < m->rows; i++) {
        for (size_t j = 0; j < m->cols; j++) {
            if (fscanf(file, "%lf", &m->data[i*m->cols + j]) != 1) {
                fclose(file);
                return -3;
            }
        }
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
 * @brief Transpose a matrix in-place.
 *
 * @param A Pointer to the matrix to transpose.
 * @return 0 on success, -1 on error.
 */
static int mat_transpose(Matrix *A) {
    if (!A || !A->data) return -1;

    size_t rows = A->rows;
    size_t cols = A->cols;
    Matrix *transposed = mat_new(cols, rows);
    if (!transposed) return -1;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            transposed->data[j * rows + i] = A->data[i * cols + j];
        }
    }

    A->rows = cols;
    A->cols = rows;
    free(A->data);
    A->data = transposed->data;

    return 0;
}


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
#define SOLVE_MAX_ITER 10000
#define SOLVE_TOL 1e-10
#define SOLVE_BRACKET_STEP 0.1
#define SOLVE_BRAC_MAX_ITER 1000


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

    // Apply permutation to b: b := P b
    double *b_copy = malloc(n * sizeof *b_copy);
    if (!b_copy) { free(perm); return LU_INVALID; }
    for (size_t i = 0; i < n; i++)
        b_copy[i] = b[perm[i]];
    for (size_t i = 0; i < n; i++)
        b[i] = b_copy[i];
    free(b_copy);

    // Forward substitution: solve L y = P b
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++)
            b[i] -= A->data[i * n + j] * b[j];
    }

    // Backward substitution: solve U x = y
    for (long i = (long)n - 1; i >= 0; i--) {
        for (size_t j = i + 1; j < n; j++)
            b[i] -= A->data[i * n + j] * b[j];
        double diag = A->data[i * n + i];
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



/**
 * @brief Solve a linear system using the Jacobi method.
 *solution stored in-place).
 * @return >=0 : number of iterations on success
 * @param A Coefficient matrix (must be square).
 * @param b Right-hand side vector (length n).
 * @param x_k Initial guess for the solution (solution stored in-place).
 * @return >=0 : number of iterations on success
 *         -1  : invalid input / memory allocation failure
 *         -2  : did not converge within max_iterations
 *         -4  : zero diagonal element encountered
 */
static int mat_jacobi_solve(const Matrix* A, const double* b, double* x_k) {
    if (!A || !b || !x_k || !A->data) return -1;
    if (A->rows != A->cols || A->rows == 0) return -1;

    size_t n = A->rows;
    double* x_next = malloc(n * sizeof(double));
    if (!x_next) return -1;

    int k = 0;
    double error;

    do {
        error = 0.0;
        for (size_t i = 0; i < n; i++) {
            double sum = 0.0;
            for (size_t j = 0; j < n; j++) {
                if (j != i) sum += A->data[i*n + j] * x_k[j];
            }
            double diag = A->data[i*n + i];
            if (fabs(diag) < 1e-14) { free(x_next); return -4; }

            double new_val = (b[i] - sum) / diag;
            error = fmax(error, fabs(new_val - x_k[i]));
            x_next[i] = new_val;
        }

        // copy back into user’s x_k
        memcpy(x_k, x_next, n * sizeof(double));

        k++;
    } while (error > SOLVE_TOL && k < SOLVE_MAX_ITER);

    free(x_next);
    return (error <= SOLVE_TOL) ? k : -2;
}



/**
 * @brief Solve a linear system using the Gauss–Seidel method with SOR.
 *
 * @param A Coefficient matrix (must be square).
 * @param b Right-hand side vector (length n).
 * @param x Initial guess for the solution (solution stored in-place).
 * @return >=0 : number of iterations on success
 *         -1  : invalid input
 *         -2  : did not converge within max_iterations
 *         -4  : zero diagonal element encountered
 */
static size_t mat_gs_solve(Matrix* A, double* b, double* x) {
    if (!A || !b || !x || !A->data) return -1;
    if (A->rows != A->cols || A->rows == 0) return -1;

    size_t n = A->rows;
    size_t k = 0;
    double error;

    double omega = 1.0;   // start as plain Gauss–Seidel
    double prev_error = INFINITY;

    do {
        error = 0.0;
        for (size_t i = 0; i < n; i++) {
            double sum = 0.0;
            for (size_t j = 0; j < n; j++) {
                if (j != i) sum += A->data[i*n + j] * x[j];
            }
            double diag = A->data[i*n + i];
            if (fabs(diag) < 1e-14) return -4;

            double raw = (b[i] - sum) / diag;
            double new_val = omega * raw + (1.0 - omega) * x[i];

            error = fmax(error, fabs(new_val - x[i]));
            x[i] = new_val;
        }

        // adaptive omega: if error stagnates, adjust relaxation
        if (error > 0.8 * prev_error) {
            omega = fmin(1.9, omega + 0.05);  // try to accelerate
        } else if (error < 0.5 * prev_error) {
            omega = fmax(1.0, omega - 0.05);  // reduce overshoot
        }
        prev_error = error;

        k++;
    } while (error > SOLVE_TOL && k < SOLVE_MAX_ITER);

    return (error <= SOLVE_TOL) ? k : -2;
}

/**
 * @brief Heuristic in-place row swapping to improve diagonal dominance of A.
 *
 * @param A Coefficient matrix (modified in-place).
 * @param b RHS vector (modified in-place).
 * @param n Matrix size.
 * @return 0 if successful, -1 if a zero diagonal could not be fixed.
 */
static int mat_diag_dom(Matrix* A, double* b) {
    if (!A || !b || !A->data) return -1;
    if (A->rows != A->cols) return -1;

    size_t n = A->rows;
    const double EPS = 1e-14;
    size_t max_passes = n;          /* tunable: number of improvement passes */

    for (size_t pass = 0; pass < max_passes; ++pass) {
        int any_swap = 0;

        for (size_t i = 0; i < n; ++i) {
            /* current row i metrics */
            double diag_i = fabs(A->data[i*n + i]);
            double offsum_i = 0.0;
            for (size_t j = 0; j < n; ++j) if (j != i) offsum_i += fabs(A->data[i*n + j]);

            /* if already (weakly) diagonally dominant for this row, skip */
            if (diag_i >= offsum_i) continue;

            /* current ratio (use -1 for unusable rows) */
            double curr_ratio = (offsum_i == 0.0) ? (diag_i > 0.0 ? INFINITY : -1.0)
                                                  : diag_i / offsum_i;

            /* search for the best row to swap into position i */
            size_t best_row = n;   /* sentinel: n means not found */
            double best_ratio = curr_ratio;

            for (size_t r = 0; r < n; ++r) {
                if (r == i) continue;

                /* if we swap row r into i, the candidate diagonal would be A[r,i] */
                double diag_r = fabs(A->data[r*n + i]);
                double offsum_r = 0.0;
                for (size_t j = 0; j < n; ++j) if (j != i) offsum_r += fabs(A->data[r*n + j]);

                double ratio;
                if (offsum_r == 0.0) {
                    ratio = (diag_r > 0.0) ? INFINITY : -1.0;
                } else {
                    ratio = diag_r / offsum_r;
                }

                /* prefer strictly better ratio (small tolerance to avoid thrash) */
                if (ratio > best_ratio + 1e-15) {
                    best_ratio = ratio;
                    best_row = r;
                }
            }

            if (best_row < n) {
                /* perform row swap i <-> best_row (in-place), and swap b */
                for (size_t c = 0; c < n; ++c) {
                    double tmp = A->data[i*n + c];
                    A->data[i*n + c] = A->data[best_row*n + c];
                    A->data[best_row*n + c] = tmp;
                }
                double tmpb = b[i];
                b[i] = b[best_row];
                b[best_row] = tmpb;
                any_swap = 1;
            } else {
                /* fallback: try a simple pivot-style swap with any row that has a bigger |A[r,i]| */
                size_t found = n;
                double max_diag = diag_i;
                for (size_t r = i + 1; r < n; ++r) {
                    double cand = fabs(A->data[r*n + i]);
                    if (cand > max_diag + 1e-15) {
                        max_diag = cand;
                        found = r;
                    }
                }
                if (found < n) {
                    for (size_t c = 0; c < n; ++c) {
                        double tmp = A->data[i*n + c];
                        A->data[i*n + c] = A->data[found*n + c];
                        A->data[found*n + c] = tmp;
                    }
                    double tmpb = b[i];
                    b[i] = b[found];
                    b[found] = tmpb;
                    any_swap = 1;
                }
            }
        } /* end loop rows i */

        if (!any_swap) break; /* stable: no improvement this pass */
    } /* end passes */

    /* Final attempt: if any diagonal is tiny/zero, try any global row that gives nonzero in this column */
    for (size_t i = 0; i < n; ++i) {
        if (fabs(A->data[i*n + i]) < EPS) {
            size_t found = n;
            for (size_t r = 0; r < n; ++r) {
                if (r == i) continue;
                if (fabs(A->data[r*n + i]) > EPS) { found = r; break; }
            }
            if (found < n) {
                for (size_t c = 0; c < n; ++c) {
                    double tmp = A->data[i*n + c];
                    A->data[i*n + c] = A->data[found*n + c];
                    A->data[found*n + c] = tmp;
                }
                double tmpb = b[i];
                b[i] = b[found];
                b[found] = tmpb;
            }
            /* if none found, we leave the zero diagonal — caller must handle (GS/SOR will fail on zero diag) */
        }
    }
    return 0;
}

/**
 * @section Section 4: Root finding algorithms for single-variable functions
 * @brief Find roots of functions using various methods.
 *
 * @details This section provides implementations of several root-finding algorithms,
 * including bisection, Newton's method, and secant method for single-variable functions. Also
 * includes Laguerre method for polynomial roots both real and complex.
 * The functions handle input validation and return codes for success or various failure modes.
 */


 /**
  * @brief Find a bracket [a, b] such that f(a) * f(b) < 0.
  * 
  * @param f Function pointer f(x).
  * @param a Pointer to the left endpoint of the interval.
  * @param b Pointer to the right endpoint of the interval.
  *
  * @return 0 on success, -1 on invalid input, -2 if no interval found.
  */
int get_bracket(double (*f)(double), double *a, double *b) {
    if (!f || !a || !b || *a >= *b) return -1;

    double fa = f(*a);
    double fb = f(*b);

    // Check if endpoints are already roots
    if (fa == 0.0) { *b = *a; return 0; }
    if (fb == 0.0) { *a = *b; return 0; }

    int iter = 0;
    while (fa * fb > 0 && iter < SOLVE_BRAC_MAX_ITER) {
        // Move the endpoint with smaller |f(x)| towards reducing |f(x)|
        if (fabs(fa) < fabs(fb)) {
            *a -= SOLVE_BRACKET_STEP*(*b -*a) ;
            fa = f(*a);
        } else {
            *b += SOLVE_BRACKET_STEP*(*b -*a) ;
            fb = f(*b);
        }
        iter++;
    }

    if (fa * fb > 0) return -2;  // Could not bracket root
    return iter;                     // Success
}


/**
 * @brief Root-finding with the Bisection method.
 *
 * @details
 * Approximates a root of f(x) in [a, b] by repeatedly halving the interval.
 * Stops when |f(c)| < SOLVE_TOL or |b - a| < SOLVE_TOL.
 *
 * @param f Function pointer f(x).
 * @param a Pointer to the left endpoint of the interval.
 * @param b Pointer to the right endpoint of the interval.
 * 
 * @return Iteration count, or:
 *   -1 invalid input, -2 no sign change, -3 max iterations hit.
 */
int root_bisection(double (*f)(double), double *a, double *b) 
{ if (!f || !a || !b || *a >= *b) return -1; 
    double fa = f(*a); 
    double fb = f(*b); 
    if (fa * fb > 0) return -2; // No root in interval
    double c, fc;
    for (int k = 0; k < SOLVE_MAX_ITER; k++) {
        c = 0.5 * (*a + *b); // midpoint
        fc = f(c);
        if (fabs(fc) < SOLVE_TOL || fabs(*b - *a) < SOLVE_TOL) {
            *a = *b = c; // collapse interval to root
            return k;
        }
        if (fa * fc < 0) {
            *b = c; fb = fc;
        } else {
            *a = c; fa = fc;
        }
    }
    return -3; // Max iterations reached
}



/**
 * @brief Root-finding with Regula Falsi + Anderson acceleration.
 *
 * @details
 * Approximates a root of f(x) in [a, b] using:
 *
 *     c = b - f(b) * (b - a) / (f(b) - f(a))
 *
 * Illinois and anderson Hybrid step: halve f(a) or f(b).
 *
 * Stops when |f(c)| < SOLVE_TOL or |b - a| < SOLVE_TOL.
 *
 * @param f Function pointer f(x).
 * @param a Pointer to the left endpoint of the interval.
 * @param b Pointer to the right endpoint of the interval.
 * 
 * @return Iteration count, or:
 *   -1 invalid input, -2 no sign change, -3 max iterations hit.
 */
int root_rf(double (*f)(double), double *a, double *b) {
    if (!f || !a || !b || *a >= *b) return -1;

    double fa = f(*a);
    double fb = f(*b);
    if (fa * fb > 0) return -2;

    double c = *a, fc;
    int side = 0; // which side updated last

    for (int k = 0; k < SOLVE_MAX_ITER; k++) {
        // regula falsi interpolation
        c = (*a * fb - *b * fa) / (fb - fa);
        fc = f(c);

        if (fabs(fc) < SOLVE_TOL || fabs(*b - *a) < SOLVE_TOL) {
            *a = *b = c;
            return k;
        }

        if (fc * fb > 0) {
            // replace b with c
            *b = c; fb = fc;
            if (side == -1) {
                // consecutive update on b → scale fa (Anderson style)
                double mprime = 1.0 - fc / fa;
                double m = (mprime > 0.0) ? mprime : 0.5;
                fa *= m;
            }
            side = -1;
        } else if (fa * fc > 0) {
            // replace a with c
            *a = c; fa = fc;
            if (side == +1) {
                // consecutive update on a → scale fb (Anderson style)
                double mprime = 1.0 - fc / fb;
                double m = (mprime > 0.0) ? mprime : 0.5;
                fb *= m;
            }
            side = +1;
        } else {
            *a = *b = c;
            return k;
        }
    }

    return -3;
}


/**
 * @brief Root finding using fixed point method (precurssor to Newton Raphson).
 * 
 * @param g Function pointer g(x) for fixed point iteration x = g(x).
 * @param x Pointer to initial guess (input), overwritten with root (output).
 * @return Iteration count on success,
 *         -1 if invalid input,
 *         -3 if not converged.
 */
static int root_fixed_point(double (*g)(double), double *x) {
    if (!g || !x) return -1;

    double xk = *x;
    double xk1;

    for (int k = 0; k < SOLVE_MAX_ITER; k++) {
        xk1 = g(xk);

        // convergence check: relative step
        if (fabs(xk1 - xk) < SOLVE_TOL * (1.0 + fabs(xk1))) {
            *x = xk1;
            return k; // converged
        }

        xk = xk1;
    }

    return -3; // failed to converge within max iterations
}



/**
 * @brief Root finding using Newton-Raphson method.
 *
 * @details
 * Iterative formula:
 *     x_{k+1} = x_k - f(x_k) / f'(x_k)
 *
 * @param f   Function pointer f(x).
 * @param df  Derivative function pointer f'(x).
 * @param x   Pointer to initial guess (input), overwritten with root (output).
 * @return Iteration count on success,
 *         -1 if invalid input,
 *         -2 if derivative too small,
 *         -3 if not converged.
 */
static int root_newton(double (*f)(double), double (*df)(double), double *x) {

    if (!f || !df || !x) return -1;

    double fx = f(*x);

    // iteration loop
    for (int k = 0; k < SOLVE_MAX_ITER; k++) {
        double dfx = df(*x);             
        if (fabs(dfx) < 1e-14) return -2; 

        double step = fx / dfx;          
        *x -= step;                       
        fx = f(*x);                       

        // convergence check: relative step or residual
        if (fabs(step) < SOLVE_TOL * (1.0 + fabs(*x)) || fabs(fx) < SOLVE_TOL)
            return k;                    
    }

    return -3; // failed to converge within max iterations
}

/**
 * @brief Root finding using the Secant method.
 *
 * @details
 * Iterative method that approximates the root of f(x) without using derivatives:
 * \f[
 *     x_{k+1} = x_k - f(x_k)\frac{x_k - x_{k-1}}{f(x_k) - f(x_{k-1})}.
 * \f]
 *
 * @param f   Function pointer f(x).
 * @param x0  Pointer to first initial guess (updated during iterations).
 * @param x1  Pointer to second initial guess; holds the final root estimate.
 *
 * @return
 *  - >= 0 : Iteration count (converged).
 *  - -1   : Invalid arguments.
 *  - -2   : Division by zero (f(x1) ≈ f(x0)).
 *  - -3   : No convergence within SOLVE_MAX_ITER.
 */
static int root_secant(double (*f)(double), double *x0, double *x1) {
    if (!f || !x0 || !x1) return -1;

    double f0 = f(*x0);
    double f1 = f(*x1);

    for (int k = 0; k < SOLVE_MAX_ITER; k++) {
        if (fabs(f1 - f0) < 1e-14) return -2; // prevent div/0

        double x2 = *x1 - f1 * (*x1 - *x0) / (f1 - f0);
        double f2 = f(x2);

        // convergence check: relative step or residual
        if (fabs(x2 - *x1) < SOLVE_TOL * (1.0 + fabs(x2)) || fabs(f2) < SOLVE_TOL) {
            *x1 = x2;
            return k; // converged
        }

        // shift for next iteration
        *x0 = *x1; f0 = f1;
        *x1 = x2; f1 = f2;
    }

    return -3; // failed to converge within max iterations
}




/**
 * @brief Find a root of a polynomial using Laguerre's method.
 *
 * @details
 * Iterates starting from *guess until convergence. Works for real or
 * complex SCALAR (depending on USE_COMPLEX).
 *
 * Formula:
 *   x_{k+1} = x_k - n / (G ± sqrt((n-1)(nH - G^2)))
 *
 * where
 *   G = p'(x)/p(x),   H = G^2 - p''(x)/p(x).
 *
 * @param P Polynomial pointer.
 * @param guess Initial guess (input), overwritten with root (output).
 * @return Iteration count if converged,
 *         -1 invalid input,
 *         -2 max iterations exceeded.
 */
int poly_root_laguerre(const Polynomial *P, SCALAR *guess) {
    if (!P || !guess || P->deg <= 0) return -1;

    SCALAR x = *guess;

    for (int k = 0; k < SOLVE_MAX_ITER; k++) {
        SCALAR px  = poly_eval(P, x);
        if (SCALAR_ABS(px) < SOLVE_TOL) {
            *guess = x;
            return k;  // converged
        }

        Polynomial *P1 = poly_derivative(P);
        Polynomial *P2 = poly_derivative(P1);

        SCALAR p1x = poly_eval(P1, x);
        SCALAR p2x = poly_eval(P2, x);

        poly_free(P1);
        poly_free(P2);

        int n = P->deg;
        SCALAR G = SCALAR_DIV(p1x, px);
        SCALAR H = SCALAR_SUB(SCALAR_MUL(G, G),
                              SCALAR_DIV(p2x, px));

        // sqrt((n-1)(nH - G^2))
        SCALAR tmp1 = SCALAR_SUB(SCALAR_MUL(CAST_SCALAR(n), H),
                                 SCALAR_MUL(G, G));
        SCALAR tmp2 = SCALAR_MUL(CAST_SCALAR(n - 1), tmp1);

        SCALAR sq = SCALAR_SQRT(tmp2);

        // choose denominator with larger magnitude
        SCALAR denom1 = SCALAR_ADD(G, sq);
        SCALAR denom2 = SCALAR_SUB(G, sq);
        SCALAR denom = (SCALAR_ABS(denom1) > SCALAR_ABS(denom2)) ? denom1 : denom2;

        if (SCALAR_ABS(denom) < 1e-14) return -2;  // prevent div/0

        SCALAR dx = SCALAR_DIV(CAST_SCALAR(n), denom);
        x = SCALAR_SUB(x, dx);

        if (SCALAR_ABS(dx) < SOLVE_TOL) {
            *guess = x;
            return k;  // converged
        }
    }

    return -2; // max iterations exceeded
}





/**
 * @section Section 5: Multivariable 
 */



 /**
 * @brief Multivariable Newton's method for solving F(x) = 0.
 *
 * @details
 * At each iteration:
 *   1. Evaluate F(x) and J(x)
 *   2. Solve J(x) Δx = -F(x) using LU decomposition
 *   3. Update x := x + Δx
 *
 * Convergence check: ||Δx|| < tol*(1+||x||) or ||F(x)|| < tol.
 *
 * @param f   Function pointer: f(x, n) → returns F(x) (length n, heap-allocated).
 * @param J   Function pointer: J(x, n) → returns Jacobian matrix (n×n).
 * @param x   Input: initial guess. Output: approximate root (length n).
 * @param n   Dimension of the system.
 *
 * @return
 *  - >= 0 : Iteration count (converged).
 *  - -1   : Invalid arguments or allocation failure.
 *  - -2   : Singular Jacobian (LU failed).
 *  - -3   : Did not converge within SOLVE_MAX_ITER.
 */
int root_multi_newton(
        double* (*f)(const double *x, int n),
        Matrix* (*J)(const double *x, int n),
        double *x,
        int n)
{
    if (!f || !J || !x || n <= 0) return -1;

    double *dx = malloc(n * sizeof(double));
    if (!dx) return -1;

    for (int k = 0; k < SOLVE_MAX_ITER; k++) {
        // Evaluate F(x)
        double *fx = f(x, n);
        if (!fx) { free(dx); return -1; }

        // Residual norm check
        double normF = 0.0;
        for (int i = 0; i < n; i++) normF += fx[i]*fx[i];
        normF = sqrt(normF);
        if (normF < SOLVE_TOL) {
            free(fx); free(dx);
            return k;
        }

        // Build RHS = -F(x)
        for (int i = 0; i < n; i++) dx[i] = -fx[i];

        // Get Jacobian
        Matrix *Jx = J(x, n);
        free(fx);
        if (!Jx) { free(dx); return -1; }

        // Solve J dx = -F(x)
        int status = mat_lu_solve(Jx, dx);
        mat_free(Jx);
        if (status != LU_SUCCESS) {
            free(dx);
            return -2;
        }

        // Update x := x + dx
        double normStep = 0.0, normX = 0.0;
        for (int i = 0; i < n; i++) {
            x[i] += dx[i];
            normStep += dx[i]*dx[i];
            normX += x[i]*x[i];
        }
        normStep = sqrt(normStep);
        normX = sqrt(normX);

        // Convergence by step size
        if (normStep < SOLVE_TOL * (1.0 + normX)) {
            free(dx);
            return k;
        }
    }

    free(dx);
    return -3; // no convergence
}






#endif  // _NUMMETH_H_
