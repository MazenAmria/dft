//
// Created by mazen on 01/02/2021.
//

#ifndef DFT_FFT_IFFT_H
#define DFT_FFT_IFFT_H
#include "Transform.h"

class IFFT : public Transform {
public:
    IFFT (real_vec data, double fs) {
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

    IFFT (complex_vec data, double fs) {
        this->fs = fs;
        try {
            init(data);
        } catch (std::exception& e) {
            throw e;
        }
    }

    void static ifft(complex_vec& data) {
        Transform::fft(data, true);
        Transform::normalize(data);
    }

private:
    void init (complex_vec& data) {
        if (data.empty()) {
            throw std::invalid_argument("Input vector cannot be empty");
        }
        this->X = data;
        this->x = this->X;
        IFFT::ifft(this->x);
    }
};


#endif //DFT_FFT_IFFT_H
