//
// Created by Mazen Amria on 01/02/2021.
//

#ifndef DFT_FFT_MATRIX_H
#define DFT_FFT_MATRIX_H
#include <vector>
#include <complex>
#include <stdexcept>


class Matrix {
    typedef std::vector<std::vector<std::complex<double>>> matrix;
    typedef std::vector<std::complex<double>> row;
    typedef std::vector<std::complex<double>> vec;
    typedef unsigned int dim;

public:
    matrix mat;
    dim shape[2]{};

    explicit Matrix (dim a) {
        try {
            init (a, a);
        } catch(std::exception& e) {
            throw e;
        }
    }

    Matrix (dim a, dim b) {
        try {
            init (a, b);
        } catch(std::exception& e) {
            throw e;
        }
    }

    static Matrix from_vector (vec data) {
        try {
            Matrix x(data.size(), 1);
            for (unsigned int i = 0; i < data.size(); i++) {
                x.mat[i][0] = data[i];
            }
            return x;
        } catch (std::exception& e) {
            throw e;
        }
    }

    vec to_vector () {
        vec data(this->shape[0] * this->shape[1]);
        for (unsigned int i = 0; i < this->shape[0]; i++) {
            for (unsigned int j = 0; j < this->shape[1]; j++) {
                data[i * this->shape[1] + j] = this->mat[i][j];
            }
        }
        return data;
    }

    Matrix add (Matrix B) {
        if (
                this->shape[0] != B.shape[0] ||
                this->shape[1] != B.shape[1]
        ) {
            throw std::invalid_argument("Addition isn't defined for Matrices with different dimensions");
        }
        Matrix result(this->shape[0], this->shape[1]);
        for (int i = 0; i < result.shape[0]; i++) {
            for (int j = 0; j < result.shape[1]; j++) {
                result.mat[i][j] = this->mat[i][j] + B.mat[i][j];
            }
        }
        return result;
    }

    Matrix subtract (Matrix B) {
        if (
                this->shape[0] != B.shape[0] ||
                this->shape[1] != B.shape[1]
                ) {
            throw std::invalid_argument("Subtraction isn't defined for Matrices with different dimensions");
        }
        Matrix result(this->shape[0], this->shape[1]);
        for (int i = 0; i < result.shape[0]; i++) {
            for (int j = 0; j < result.shape[1]; j++) {
                result.mat[i][j] = this->mat[i][j] - B.mat[i][j];
            }
        }
        return result;
    }

    Matrix multiply (Matrix B) {
        if (this->shape[1] != B.shape[0]) {
            throw std::invalid_argument("2nd Dimension of left matrix should equal the 1st dimension of the right");
        }
        Matrix result(this->shape[0], B.shape[1]);
        for (int i = 0; i < result.shape[0]; i++) {
            for (int j = 0; j < result.shape[1]; j++) {
                result.mat[i][j] = 0;
                for(int k = 0; k < B.shape[0]; k++) {
                    result.mat[i][j] += this->mat[i][k] * B.mat[k][j];
                }
            }
        }
        return result;
    }

    Matrix transpose () {
        Matrix result(this->shape[1], this->shape[0]);
        for (int i = 0; i < result.shape[0]; i++) {
            for (int j = 0; j < result.shape[1]; j++) {
                result.mat[i][j] = this->mat[j][i];
            }
        }
        return result;
    }

    int determinant () {
        if (this->shape[0] != this->shape[1]) {
            throw std::invalid_argument("Determinant isn't defined for Non-Square Matrix");
        }
        // TODO: implement Matrix.determinant () {...}
        return 0;
    }

    Matrix inverse () {
        try {
            if (this->determinant() == 0) {
                throw std::invalid_argument("Inverse isn't defined for Singular Matrix");
            }
            // TODO: implement Matrix.inverse () {...}
            return Matrix(1);
        } catch (std::exception& e) {
            throw e;
        }
    }

private:
    void init (dim a, dim b) {
        if (a == 0 || b == 0) {
            throw std::invalid_argument("Matrix cannot have zero dimensions");
        }
        this->mat = matrix(a, row(b));
        this->shape[0] = a;
        this->shape[1] = b;
    }
};


#endif //DFT_FFT_MATRIX_H
