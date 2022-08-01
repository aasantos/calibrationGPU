#!/bin/bash
nvcc -ptx -std=c++11 -O2 -arch=sm_86 -I/home/antonio/Develop/CUDA/include bates.cu
