#pragma once
#include <vector>
#include <Eigen/Core>
#include <utility>

template <typename Scalar, int Rows, int Cols>
using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;

template <typename Scalar, int Rows>
using Vector = Eigen::Matrix<Scalar, Rows, 1>;

// Diagonal Matrix
template <typename Scalar, int Dimisions>
class DiagonalMatrix
{
public:
    DiagonalMatrix() = default;
    DiagonalMatrix(const DiagonalMatrix& rhs) : data_(rhs.data_) {};
    DiagonalMatrix(DiagonalMatrix&& rhs) : data_(std::move(rhs.data_)) {};

    operator Matrix<Scalar, Dimisions, Dimisions>()
    {
        Matrix<Scalar, Dimisions, Dimisions> res;
        DiagonalMatrix<Scalar, Dimisions>& self = *this;
        res.setZero();
        for (int i = 0; i < Dimisions; ++i)
        {
            res(i, i) = self(i);
        }
        return res;
    }

    // Diagonl Matrix operator overload
    DiagonalMatrix operator+(const DiagonalMatrix& rhs) const
    {
        DiagonalMatrix res;
        const DiagonalMatrix& self = *this;
        for (int i = 0; i < Dimisions; ++i)
        {
            res(i) = self(i) + rhs(i);
        }
        return res;
    }
    DiagonalMatrix& operator+=(const DiagonalMatrix& rhs)
    {
        DiagonalMatrix& self = *this;
        for (int i = 0; i < Dimisions; ++i)
        {
            self(i) += rhs(i);
        }
        return self;
    }
    DiagonalMatrix operator-() const
    {
        DiagonalMatrix res;
        const DiagonalMatrix& self = *this;
        for (int i = 0; i < Dimisions; ++i)
        {
            res(i) = -self(i);
        }
        return res;
    }
    DiagonalMatrix operator-(const DiagonalMatrix& rhs) const
    {
        DiagonalMatrix res;
        const DiagonalMatrix& self = *this;
        for (int i = 0; i < Dimisions; ++i)
        {
            res(i) = self(i) - rhs(i);
        }
        return res;
    }
    DiagonalMatrix& operator-=(const DiagonalMatrix& rhs)
    {
        DiagonalMatrix& self = *this;
        for (int i = 0; i < Dimisions; ++i)
        {
            self(i) -= rhs(i);
        }
        return self;
    }
    DiagonalMatrix& operator=(const DiagonalMatrix& rhs)
    {
        // self assignment check
        if (this != &rhs)
        {
            data_ = rhs.data_;
        }
        return *this;
    }
    DiagonalMatrix operator*(const DiagonalMatrix& rhs) const
    {
        DiagonalMatrix res;
        for (int i = 0; i < Dimisions; ++i)
        {
            res.data_(i) = data_(i) * rhs.data_(i);
        }
        return res;
    }
    DiagonalMatrix& operator*=(const DiagonalMatrix& rhs)
    {
        for (int i = 0; i < Dimisions; ++i)
        {
            data_(i) = data_(i) * rhs.data_(i);
        }
        return *this;
    }

    Scalar operator()(const int index) const
    {
        assert(index >= 0 && index < Dimisions);
        return data_(index);
    }
    Scalar& operator()(const int index)
    {
        assert(index >= 0 && index < Dimisions);
        return data_(index);
    }

private:
    Vector<Scalar, Dimisions> data_;
};

// Diagnal Matrix operator overload
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator*(const DiagonalMatrix<Scalar, Dimisions>& lhs, const Matrix<Scalar, Dimisions, Dimisions>& rhs)
{
    Matrix<Scalar, Dimisions, Dimisions> res;
    for (int i = 0; i < Dimisions; ++i)
    {
        res.rows(i) = rhs.rows(i) * lhs(i);
    }
    return res;
}

// Diagnal Matrix operator overload
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator*(const Matrix<Scalar, Dimisions, Dimisions>& lhs, const DiagonalMatrix<Scalar, Dimisions>& rhs)
{
    Matrix<Scalar, Dimisions, Dimisions> res;
    for (int i = 0; i < Dimisions; ++i)
    {
        res.cols(i) = lhs.cols(i) * rhs(i);
    }
    return res;
}

template <typename Scalar, int Dimisions>
Vector<Scalar, Dimisions> operator*(const Matrix<Scalar, 1, Dimisions>& lhs, const DiagonalMatrix<Scalar, Dimisions>& rhs)
{
    Vector<Scalar, Dimisions> res;
    for (int i = 0; i < Dimisions; ++i)
    {
        res(i) = lhs(i) * rhs(i);
    }
}

template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator+(const DiagonalMatrix<Scalar, Dimisions>& lhs, const Matrix<Scalar, Dimisions, Dimisions>& rhs)
{
    Matrix<Scalar, Dimisions, Dimisions> res = rhs;
    for (int i = 0; i < Dimisions; ++i)
    {
        res(i, i) += lhs(i);
    }
    return res;
}
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator+(const Matrix<Scalar, Dimisions, Dimisions>& lhs, const DiagonalMatrix<Scalar, Dimisions>& rhs)
{
    return rhs + lhs;
}

template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator-(const Matrix<Scalar, Dimisions, Dimisions>& lhs, const DiagonalMatrix<Scalar, Dimisions>& rhs)
{
    Matrix<Scalar, Dimisions, Dimisions> res = lhs;
    for (int i = 0; i < Dimisions; ++i)
    {
        res(i, i) -= rhs(i);
    }
    return res;
}

template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator-(const DiagonalMatrix<Scalar, Dimisions>& lhs, const Matrix<Scalar, Dimisions, Dimisions>& rhs)
{
    Matrix<Scalar, Dimisions, Dimisions> res = -rhs;
    for (int i = 0; i < Dimisions; ++i)
    {
        res(i, i) += lhs(i);
    }
    return res;
}
