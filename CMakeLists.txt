cmake_minimum_required(VERSION 3.17)
project(DFT_FFT)

set(CMAKE_CXX_STANDARD 20)

add_executable(DFT_FFT main.cpp transform/Matrix.h transform/DFT.h transform/IDFT.h transform/FFT.h transform/IFFT.h Transform.cpp transform/Transform.h)

find_package(Python COMPONENTS Development)
target_include_directories(DFT_FFT PRIVATE ${Python_INCLUDE_DIRS})
target_link_libraries(DFT_FFT Python::Python)