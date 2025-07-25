#include "vec.cuh"

// simple 5D mathematical vectors that can be added and multiplied
__host__ __device__ RGRP::kernels::vec5d::vec5d() : data{0.0, 0.0, 0.0, 0.0, 0.0} {}
__host__ __device__ RGRP::kernels::vec5d::vec5d(real x0, real x1, real x2, real x3, real x4) : data{x0, x1, x2, x3, x4} {}

__host__ __device__ real &RGRP::kernels::vec5d::operator[](size_t i) { return data[i]; }
__host__ __device__ real &RGRP::kernels::vec5d::operator()(size_t i) { return data[i]; }
__host__ __device__ const real &RGRP::kernels::vec5d::operator[](size_t i) const { return data[i]; }
__host__ __device__ const real &RGRP::kernels::vec5d::operator()(size_t i) const { return data[i]; }

__host__ __device__ RGRP::kernels::vec5d RGRP::kernels::vec5d::operator+(const RGRP::kernels::vec5d &v)
{
    vec5d ret;
    for (auto i = 0; i < 4; i++)
    {
        ret[i] = this->data[i] + v[i];
    }
    return ret;
}

__host__ __device__ RGRP::kernels::vec5d RGRP::kernels::vec5d::operator-(const RGRP::kernels::vec5d &v)
{
    vec5d ret;
    for (auto i = 0; i < 4; i++)
    {
        ret[i] = this->data[i] - v[i];
    }
    return ret;
}

__host__ __device__ RGRP::kernels::vec5d RGRP::kernels::vec5d::operator+(const RGRP::kernels::vec5d &v) const
{
    vec5d ret;
    for (auto i = 0; i < 4; i++)
    {
        ret[i] = this->data[i] + v[i];
    }
    return ret;
}

__host__ __device__ RGRP::kernels::vec5d RGRP::kernels::vec5d::operator-(const RGRP::kernels::vec5d &v) const
{
    vec5d ret;
    for (auto i = 0; i < 4; i++)
    {
        ret[i] = this->data[i] - v[i];
    }
    return ret;
}

__host__ __device__ RGRP::kernels::vec5d RGRP::kernels::vec5d::operator*(real a)
{
    vec5d ret;
    for (auto i = 0; i < 4; i++)
    {
        ret[i] = this->data[i] * a;
    }
    return ret;
}
