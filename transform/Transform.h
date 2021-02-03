//
// Created by Mazen Amria on 02/02/2021.
//

#ifndef DFT_FFT_TRANSFORM_H
#define DFT_FFT_TRANSFORM_H
#include "Matrix.h"


class Transform {
protected:
    typedef std::vector<std::complex<double>> complex_vec;
    typedef std::vector<double> real_vec;
    typedef std::pair<std::vector<double>, std::vector<double>> plot;
    typedef std::complex<double> complex_d;

    complex_vec x;
    complex_vec X;
    double fs;

    void static dft(complex_vec& data, std::complex<double> w) {
        try {
            Matrix X_wrap = Matrix::from_vector(data);
            Matrix transformation(data.size());
            std::complex<double> w_ = 1;
            for (unsigned int i = 0; i < data.size(); i++) {
                std::complex<double> _w_ = 1;
                for (unsigned int j = 0; j < data.size(); j++) {
                    transformation.mat[i][j] = _w_;
                    transformation.mat[j][i] = transformation.mat[i][j];
                    _w_ *= w_;
                }
                w_ *= w;
            }
            Matrix x_wrap = transformation.multiply(X_wrap);
            data = x_wrap.to_vector();
        } catch (std::exception& e) {
            throw e;
        }
    }

    void static fft(complex_vec& data, bool inverse) {
        unsigned int n = data.size();
        if ((n & (n - 1)) != 0) {
            Transform::ensure_length(data);
            n = data.size();
        }
        if (n == 1) {
            return;
        }
        complex_vec even(n / 2), odd(n / 2);
        for (int i = 0; i < n / 2; i++) {
            even[i] = data[2 * i];
            odd[i] = data[2 * i + 1];
        }
        Transform::fft(even, inverse);
        Transform::fft(odd, inverse);

        complex_d w_ = 1, w = exp(complex_d(0.0,2.0 * M_PI / (double) data.size()));
        for (unsigned int i = 0; i < data.size() / 2; i++) {
            data[i] = even[i] + w_ * odd[i];
            data[i + n / 2] = even[i] - w_ * odd[i];
            w_ *= w;
        }
    }

public:
    complex_vec transformed () {
        return this->X;
    }

    complex_vec i_transformed () {
        return this->x;
    }

    real_vec i_real_part () {
        complex_vec data = this->x;
        real_vec real_data(data.size());
        for (unsigned int i = 0; i < data.size(); i++) {
            real_data[i] = std::real(data[i]);
        }
        return real_data;
    }

    real_vec i_imaginary_part () {
        complex_vec data = this->x;
        real_vec imag_data(data.size());
        for (unsigned int i = 0; i < data.size(); i++) {
            imag_data[i] = std::imag(data[i]);
        }
        return imag_data;
    }

    real_vec real_part () {
        complex_vec data = this->X;
        real_vec real_data(data.size());
        for (unsigned int i = 0; i < data.size(); i++) {
            real_data[i] = std::real(data[i]);
        }
        return real_data;
    }

    real_vec imaginary_part () {
        complex_vec data = this->X;
        real_vec imag_data(data.size());
        for (unsigned int i = 0; i < data.size(); i++) {
            imag_data[i] = std::imag(data[i]);
        }
        return imag_data;
    }

    plot amplitude_spectrum () {
        complex_vec data = this->X;
        real_vec amplitude(data.size());
        for (unsigned int i = 0; i < data.size(); i++) {
            amplitude[i] = std::abs(data[i]);
            if (amplitude[i] < 1e-6 && amplitude[i] > -1e-6)
                amplitude[i] = 0;
        }
        Transform::shift(amplitude);
        Transform::normalize(amplitude);
        return {this->shifted_frequency(), amplitude};
    }

    plot phase_spectrum () {
        complex_vec data = this->X;
        real_vec phase(data.size());
        for (unsigned int i = 0; i < data.size(); i++) {
            phase[i] = std::arg(data[i]);
        }
        Transform::shift(phase);
        Transform::unwrap(phase);
        return {this->shifted_frequency(), phase};
    }

    real_vec time () {
        unsigned int n = this->x.size();
        real_vec t(n);
        for (unsigned int i = 0; i < n; i++) {
            t[i] = i / this->fs;
        }
        return t;
    }

    real_vec frequency () {
        unsigned int n = this->X.size();
        real_vec f(n);
        for (unsigned int i = 0; i < n; i++) {
            f[i] = i * this->fs / n;
        }
        return f;
    }

    real_vec shifted_frequency () {
        unsigned int n = this->X.size(), j = 0;
        real_vec f(n);
        for (int i = (int)(n - 1) / -2; j < n; j++, i++) {
            f[j] = i * this->fs / n;
        }
        return f;
    }

    template <class T>
    void static shift(T& data) {
        unsigned int nq_rate = data.size() / 2;
        data.insert(
                data.end(),
                data.begin(),
                data.begin() + nq_rate + 1
        );
        data.erase(
                data.begin(),
                data.begin() + nq_rate + 1
        );
    }

    template <class T>
    void static normalize(T& data) {
        for (unsigned int i = 0; i < data.size(); i++) {
            data[i] /= data.size();
        }
    }

    void static unwrap(real_vec& data) {
        // TODO: check Transform::unwrap (data) {...}
        real_vec new_data(data.size());
        new_data[(data.size() - 1) / 2] = 0;
        for (unsigned int i = (data.size() - 1) / 2; i > 0; i--) {
            double d = data[i - 1] - data[i];
            d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
            new_data[i - 1] = new_data[i] + d;
        }
        for (unsigned int i = (data.size() + 1) / 2; i < data.size(); i++) {
            double d = data[i] - data[i - 1];
            d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
            new_data[i] = new_data[i - 1] + d;
        }
    }

    void static ensure_length (complex_vec& data) {
        unsigned int appropriate_length = 1;
        unsigned int base = 2;
        unsigned int n = (unsigned int) ceil(log2(data.size()));
        while (n > 0) {
            if (n & 1)
                appropriate_length *= base;
            n /= 2;
            base *= base;
        }
        complex_vec new_data(appropriate_length - data.size());
        new_data.insert(
                new_data.end(),
                data.begin(),
                data.end()
        );
        data = new_data;
    }
};


#endif //DFT_FFT_TRANSFORM_H
