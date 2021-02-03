//
// Created by mazen on 01/02/2021.
//

#ifndef DFT_FFT_FFT_H
#define DFT_FFT_FFT_H
#include "Transform.h"


class FFT : public Transform {
public:
    FFT (real_vec data, double fs) {
        complex_vec data_wrap(data.size());
        for (unsigned int i = 0; i < data.size(); i++) {
            data_wrap[i] = data[i];
        }
        this->fs = fs;
        try {
            init(data_wrap);
        } catch (std::exception& e) {
            throw e;
        }
    }

    FFT (complex_vec data, double fs) {
        this->fs = fs;
        try {
            init(data);
        } catch (std::exception& e) {
            throw e;
        }
    }

    void static fft(complex_vec& data) {
        std::complex<double> w = exp(std::complex<double>(0.0,2.0 * M_PI / (double) data.size()));
        Transform::fft(data, w);
    }

private:
    void init (complex_vec& data) {
        if (data.empty()) {
            throw std::invalid_argument("Input vector cannot be empty");
        }
        this->X = this->x;
        FFT::fft(this->X);
    }
};


#endif //DFT_FFT_FFT_H
