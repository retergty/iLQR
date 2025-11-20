#pragma once

#include "Types.hpp"
#include <Eigen/LU>
#include <Eigen/Cholesky>

namespace LinearAlgebra
{
    /**
     * Computes the U*U^T decomposition associated to the inverse of the input matrix, where U is an upper triangular
     * matrix. Note that the U*U^T decomposition is different from the Cholesky decomposition (U^T*U).
     *
     * @param [in] Am: A symmetric square positive definite matrix
     * @param [out] AmInvUmUmT: The upper-triangular matrix associated to the UUT decomposition of inv(Am) matrix.
     */
    template <typename Scalar, int Dimisions>
    void computeInverseMatrixUUT(const Matrix<Scalar, Dimisions, Dimisions> &Am, Matrix<Scalar, Dimisions, Dimisions> &AmInvUmUmT)
    {
        // Am = Lm Lm^T --> inv(Am) = inv(Lm^T) inv(Lm) where Lm^T is upper triangular
        Eigen::LLT<Matrix<Scalar, Dimisions, Dimisions>> lltOfA(Am);
        AmInvUmUmT.setIdentity(Am.rows(), Am.cols()); // for dynamic size matrices
        lltOfA.matrixU().solveInPlace(AmInvUmUmT);
    }
}