#!/bin/bash
nvcc -ptx -std=c++11 -arch=sm_86 -I/home/antonio/Develop/CUDA/include heston.cu
