#include <iostream>
#include <vector>
#include "Transform.cpp"
#include "includes/exprtk.hpp"
#include "includes/matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
    double fs, t;
    unsigned int n;
    std::string user_input;
    expr:
    std::cin.clear();
    fflush(stdin);
    std::cout << "Expression: ";
    getline(std::cin, user_input);

    exprtk::symbol_table<double> _symbol_table;
    _symbol_table.add_variable("t",t);
    _symbol_table.add_constants();

    exprtk::expression<double> _expression;
    _expression.register_symbol_table(_symbol_table);

    exprtk::parser<double> _parser;
    if (!_parser.compile(user_input, _expression)) {
        std::cerr << "Invalid Expression." << std::endl;
        goto expr;
    }

    std::cin.clear();
    fflush(stdin);
    std::cout << "Sampling frequency: ";
    std::cin >> fs;

    std::cin.clear();
    fflush(stdin);
    std::cout << "Number of samples: ";
    std::cin >> n;

    std::vector<double> x(n);
    for (unsigned int i = 0; i < n; i++) {
        t = (double) i / fs;
        x[i] = _expression.value();
    }
    try {
        FFT transform(x, fs);
        IFFT i_transform(transform.transformed(), fs);
        std::pair<std::vector<double>, std::vector<double>> data;
        while (true) {
            std::cout << "Select: " << std::endl;
            std::cout << "1- Plot the function" << std::endl;
            std::cout << "2- Plot the Amplitude Spectrum" << std::endl;
            std::cout << "3- Plot the Phase Spectrum" << std::endl;
            unsigned int choice;
            std::cin >> choice;
            switch (choice) {
                case 1:
                    plt::figure();
                    plt::plot(i_transform.time(), i_transform.i_real_part());
                    plt::show();
                    break;
                case 4:
                    plt::figure();
                    plt::plot(transform.time(), transform.i_real_part());
                    plt::show();
                    break;
                case 2:
                    plt::figure();
                    data = transform.amplitude_spectrum();
                    plt::plot(data.first, data.second);
                    plt::show();
                    break;
                case 3:
                    plt::figure();
                    data = transform.phase_spectrum();
                    plt::plot(data.first, data.second);
                    plt::show();
                    break;
                default:
                    std::cout << "Invalid Selection" << std::endl;
            }
        }
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}