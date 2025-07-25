#pragma once
#include <cuda_runtime.h>
using real = double;

namespace RGRP::kernels
{
    // simple 5D mathematical vectors that can be added and multiplied
    struct vec5d
    {
        __host__ __device__ vec5d();
        __host__ __device__ vec5d(real x0, real x1, real x2, real x3, real x4);
        __host__ __device__ real &operator[](size_t i);
        __host__ __device__ real &operator()(size_t i);
        __host__ __device__ const real &operator[](size_t i) const;
        __host__ __device__ const real &operator()(size_t i) const;
        __host__ __device__ vec5d operator+(const vec5d &v) const;
        __host__ __device__ vec5d operator-(const vec5d &v) const;
        __host__ __device__ vec5d operator+(const vec5d &v);
        __host__ __device__ vec5d operator-(const vec5d &v);
        __host__ __device__ vec5d operator*(real a);
        real data[5];
    };
}