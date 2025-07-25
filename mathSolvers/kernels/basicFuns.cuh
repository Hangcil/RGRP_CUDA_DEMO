#pragma once
#include "vec.cuh"
#include <math.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
using thrust::device_vector;
using thrust::host_vector;
using thrust::max;
using thrust::min;
using thrust::reduce;
using thrust::swap;
using RGRP::kernels::vec5d;

namespace RGRP::kernels
{
    // flux function at a cell interface
    __device__ vec5d flux(const vec5d &U, const vec5d &U_t, real k);

    // MUSCL reconstruction procedure
    __global__ void MUSCL(vec5d *U, vec5d *Us, int Nx, int Ny, real alpha, real width, vec5d *U_lr_interface);

    // compute the approximate Riemann solutions and correct the negativity
    __global__ void rel_RP_posi_fix(vec5d *U, vec5d *Us, vec5d *U_lr_interface, int Nx, int Ny, real gamma);
    
    __global__ void GRP_l(vec5d *U, vec5d *Us, vec5d *U_lr_interface, int Nx, int Ny, real gamma);

    // convert (\rho,p,u,v) to (\rho,\rho u,\rho u^2+p,\rho v)
    __device__ vec5d rhoPUVToConserVar(const vec5d &rhoPUV, real gamma);

    // convert (\rho,\rho u,\rho u^2+p,\rho v) to (\rho,p,u,v)
    __device__ vec5d conserVarToRhoPUV(const vec5d &U_, real gamma);

    // compute the maximum characteristic speed
    __global__ void maxPropagtingSpeed(vec5d *U, real *speed, int Nx, int Ny, real gamma);

    // update using FVM method
    __global__ void forward_x_dir(vec5d *U, vec5d *U_lr_interface, int Nx, int Ny, real width, real t, real gamma);

    // swap directions
    __global__ void rotate_and_flip_vertically(vec5d *input, vec5d *output, int Nx, int Ny);

    __global__ void combineResults(vec5d *U1, vec5d *U2, int Nx, int Ny);
}//end namespace RGRP::computeKernels