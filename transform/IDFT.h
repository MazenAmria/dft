//
// Created by Mazen Amria on 01/02/2021.
//

#ifndef DFT_FFT_IDFT_H
#define DFT_FFT_IDFT_H
#include "Transform.h"


class IDFT : public Transform {
public:
    IDFT (real_vec data, double fs) {
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

    IDFT (complex_vec data, double fs) {
        this->fs = fs;
        try {
            init(data);
        } catch (std::exception& e) {
            throw e;
        }
    }

    void static idft(complex_vec& data) {
        try {
            Transform::dft(data, exp(std::complex<double>(0.0,2.0 * M_PI / (double) data.size())));
        } catch (std::exception& e) {
            throw e;
        }
    }

private:
    void init (complex_vec& data) {
        if (data.empty()) {
            throw std::invalid_argument("Input vector cannot be empty");
        }
        this->X = data;
        this->x = data;
        try {
            IDFT::idft(this->x);
        } catch (std::exception& e) {
            throw e;
        }
    }
};


#endif //DFT_FFT_IDFT_H
