#pragma once
#include <vector>
#include <cmath>
#include <mkl.h>
#include <stdexcept>
#include <iostream>
#include <omp.h>
#include "globals/mpi_vars.hpp"
#include <iostream>

template <typename XContainer, typename YContainer>
inline double integrate(const XContainer& x, const YContainer& y) {
    size_t Nx = x.size();
    size_t Ny = y.size();

    if (Ny == Nx) {
        // Trapezoidal rule
        double result = 0.0;
        for (size_t i = 0; i < Nx - 1; ++i) {
            double dx = x[i+1] - x[i];
            result += 0.5 * (y[i] + y[i+1]) * dx;
        }
        return result;
    } else if (Ny == Nx - 1) {
        // Left Riemann sum
        double result = 0.0;
        for (size_t i = 0; i < Ny; ++i) {
            double dx = x[i+1] - x[i];
            result += y[i] * dx;
        }
        return result;
    } else {
        throw std::invalid_argument("integrate(): y.size() must be equal to x.size() or x.size() - 1.");
    }
}

inline std::vector<double> solveLeastSquaresQR_MKL(const std::vector<double>& A_input,
    const std::vector<double>& b_input,
    int m, int n) {
    // Copy A and b to work arrays since LAPACK routines modify them in-place
    std::vector<double> A = A_input;          // Size m × n
    std::vector<double> b = b_input;          // Size m
    // if (A.size() != size_t(m * n) || b.size() != size_t(m)) {
    // throw std::invalid_argument("Matrix/vector size mismatch.");
    // }

    // LAPACKE_dgels expects leading dimension >= max(1, m)
    lapack_int lda = m;
    lapack_int ldb = std::max(m, n);  // b will be padded to length ldb

    // Resize b to match ldb
    b.resize(ldb, 0.0);

    lapack_int info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', m, n, 1,
    A.data(), lda,
    b.data(), ldb);

    if (info != 0) {
    throw std::runtime_error("LAPACKE_dgels failed with error code " + std::to_string(info));
    }

    // Solution x is in the first n entries of b
    return std::vector<double>(b.begin(), b.begin() + n);
}

inline std::vector<double> solveNormalEquationManual(
    const std::vector<double>& A,
    const std::vector<double>& b,
    int m, int n)
{
    // if (A.size() != size_t(m * n) || b.size() != size_t(m)) {
    //     throw std::invalid_argument("Matrix/vector size mismatch.");
    // }

    // Step 1: Compute AtA = Aᵗ * A (n×n)
    std::vector<double> AtA(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < m; ++k) {
                sum += A[k + i * m] * A[k + j * m];  // column-major access
            }
            AtA[i + j * n] = sum;
            if (i != j) AtA[j + i * n] = sum;  // Symmetric fill
        }
    }

    // Step 2: Compute Atb = Aᵗ * b (n)
    std::vector<double> Atb(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = 0; k < m; ++k) {
            sum += A[k + i * m] * b[k];
        }
        Atb[i] = sum;
    }

    // Step 3: Cholesky decomposition: AtA = L * Lᵗ
    std::vector<double> L(n * n, 0.0);  // Lower triangular matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = AtA[i + j * n];
            for (int k = 0; k < j; ++k)
                sum -= L[i + k * n] * L[j + k * n];

            if (i == j) {
                if (sum <= 0.0) throw std::runtime_error("Matrix is not positive definite.");
                L[i + j * n] = std::sqrt(sum);
            } else {
                L[i + j * n] = sum / L[j + j * n];
            }
        }
    }

    // Step 4: Solve L y = Atb (forward substitution)
    std::vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = Atb[i];
        for (int j = 0; j < i; ++j) {
            sum -= L[i + j * n] * y[j];
        }
        y[i] = sum / L[i + i * n];
    }

    // Step 5: Solve Lᵗ x = y (back substitution)
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = y[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= L[j + i * n] * x[j];
        }
        x[i] = sum / L[i + i * n];
    }

    return x;
}

inline std::vector<double> solveNormalEquationGauss(
    const std::vector<double>& A,
    const std::vector<double>& b,
    int m, int n)
{
    // if (A.size() != size_t(m * n) || b.size() != size_t(m)) {
    //     throw std::invalid_argument("Matrix/vector size mismatch.");
    // }

    // Step 1: Compute AtA = Aᵗ * A (n×n) and Atb = Aᵗ * b (n)
    std::vector<double> AtA(n * n, 0.0);
    std::vector<double> Atb(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < m; ++k) {
                sum += A[k + i * m] * A[k + j * m];
            }
            AtA[i + j * n] = sum;
            AtA[j + i * n] = sum;  // symmetric
        }

        double sum_b = 0.0;
        for (int k = 0; k < m; ++k) {
            sum_b += A[k + i * m] * b[k];
        }
        
        Atb[i] = sum_b;
    }
    // Step 2: Gaussian elimination on AtA and Atb
    for (int k = 0; k < n - 1; ++k) {
        double pivot = AtA[k + k * n];
        
        if (pivot == 0.0) {
            throw std::runtime_error("Zero pivot encountered in Gaussian elimination.");
        }

        for (int i = k + 1; i < n; ++i) {
            double factor = AtA[i + k * n] / pivot;
            
            AtA[i + k * n] = factor;  // store L

            for (int j = k + 1; j < n; ++j) {
                AtA[i + j * n] -= factor * AtA[k + j * n];
            }
        }
    }
    
    // Step 3: Forward substitution (solve L * y = Atb)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            Atb[i] -= AtA[i + j * n] * Atb[j];
        }
    }
    // Step 4: Backward substitution (solve U * x = y)
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            Atb[i] -= AtA[i + j * n] * Atb[j];
        }

        double diag = AtA[i + i * n];
        if (diag == 0.0) {
            throw std::runtime_error("Zero diagonal encountered in back substitution.");
        }
        Atb[i] /= diag;
    }

    // Atb now contains the solution vector x
    return Atb;
}


