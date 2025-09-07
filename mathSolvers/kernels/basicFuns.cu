#include "basicFuns.cuh"

// flux function at a cell interface
__device__ vec5d RGRP::kernels::flux(const vec5d &U, const vec5d &U_t, real k)
{
    real rho = U(0);
    real pi = U(1);
    real u = U(2);
    real v = U(3);
    real e = U(4);

    real rho_t = U_t(0);
    real pi_t = U_t(1);
    real u_t = U_t(2);
    real v_t = U_t(3);
    real e_t = U_t(4);

    return {
        rho * u + k / 2.0 * (rho * u_t + u * rho_t),
        rho * u * u + pi + k / 2.0 * (rho_t * u * u + 2.0 * rho * u * u_t + pi_t),
        rho * u * (e + u * u / 2.0 + v * v / 2.0) + pi * u + k / 2.0 * (rho_t * u * e + rho * u_t * e + rho * u * e_t + rho_t * u * u * u / 2.0 + 1.5 * rho * u * u * u_t + rho_t * u * v * v / 2.0 + rho * u * v * v_t + rho * u_t * v * v / 2.0 + pi_t * u + pi * u_t),
        rho * u * v + k / 2.0 * (rho * u_t * v + u * rho_t * v + u * rho * v_t),
        0.0};
}

// MUSCL reconstruction procedure
__global__ void RGRP::kernels::MUSCL(vec5d *U, vec5d *Us, int Nx, int Ny, real alpha, real width, vec5d *U_lr_interface)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index
    real h = width / real(Nx);

    if (i < Nx - 1 && j < Ny && i > 0)
    {
        int idxL = j * Nx + i - 1; // left cell
        int idx = j * Nx + i;      // cell i
        int idxR = j * Nx + i + 1; // Right cell

        vec5d Ul_s = (U[idxR] - U[idxL]) * 0.5 * (1.0 / h);
        vec5d U_s = (U[idx] - U[idxL]) * (1.0 / h) * alpha;
        vec5d Ur_s = (U[idxR] - U[idx]) * (1.0 / h) * alpha;
        for (auto i = 0; i <= 3; i++)
        {
            real sgn13 = Ul_s(i) * Ur_s(i), sgn23 = U_s(i) * Ur_s(i);
            if (sgn13 > 0.0 && sgn23 > 0.0)
            {
                real sgn = 1.0;
                if (Ul_s(i) <= 0.0)
                {
                    sgn = -1.0;
                }
                Us[idx][i] = sgn * min(abs(Ul_s(i)), min(abs(U_s(i)), abs(Ur_s(i))));
            }
            else
            {
                Us[idx][i] = 0.0;
            }
        }
        U_lr_interface[2 * idx + 1] = U[idx] - Us[idx] * h * 0.5;
        U_lr_interface[2 * idx + 2] = U[idx] + Us[idx] * h * 0.5;
    }
}

// compute the approximate Riemann solutions and correct the negativity
__global__ void RGRP::kernels::rel_RP_posi_fix(vec5d *U, vec5d *Us, vec5d *U_lr_interface, int Nx, int Ny, real gamma)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index

    if (i < Nx - 1 && j < Ny && i > 1)
    {
        int idxL = j * Nx + i - 1; // left cell
        int idx = j * Nx + i;      // cell i(left interface)

        // The implementation miantains proprietary untill the date of submission of my paper.
    }
}

// convert (\rho,p,u,v) to (\rho,\rho u,\rho u^2+p,\rho v)
__device__ vec5d RGRP::kernels::rhoPUVToConserVar(const vec5d &rhoPUV, real gamma)
{
    return {rhoPUV(0),
            rhoPUV(0) * rhoPUV(2),
            rhoPUV(1) / (gamma - 1.0) + rhoPUV(0) * (rhoPUV(2) * rhoPUV(2) + rhoPUV(3) * rhoPUV(3)) / 2.0,
            rhoPUV(0) * rhoPUV(3),
            0.0};
}

// convert (\rho,\rho u,\rho u^2+p,\rho v) to (\rho,p,u,v)
__device__ vec5d RGRP::kernels::conserVarToRhoPUV(const vec5d &U_, real gamma)
{
    return {U_(0),
            (gamma - 1.0) * (U_(2) - U_(1) * U_(1) / U_(0) / 2.0 - U_(3) * U_(3) / U_(0) / 2.0),
            U_(1) / U_(0),
            U_(3) / U_(0),
            0.0};
}

// compute the maximum characteristic speed
__global__ void RGRP::kernels::maxPropagtingSpeed(vec5d *U, real *speed, int Nx, int Ny, real gamma)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index

    if (i < Nx && j < Ny)
    {
        int idx = j * Nx + i;
        real velo = sqrt(U[idx](2) * U[idx](2) + U[idx](3) * U[idx](3)), rho = U[idx](0), p = U[idx](1), c = sqrt(gamma * p / rho);
        speed[idx] = abs(velo) + c;
    }
}

// update using FVM method
__global__ void RGRP::kernels::forward_x_dir(vec5d *U, vec5d *U_lr_interface, int Nx, int Ny, real width, real t, real gamma)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index
    real h = width / real(Nx);

    if (i < Nx - 2 && j < Ny && i > 1)
    {
        int idx = j * Nx + i; // cell i
        vec5d F_L = flux(U_lr_interface[2 * idx], U_lr_interface[2 * idx + 1], t);
        vec5d F_R = flux(U_lr_interface[2 * idx + 2], U_lr_interface[2 * idx + 3], t);
        vec5d U_ = rhoPUVToConserVar(U[idx], gamma);
        real lambda = t / h;
        U_ = U_ - (F_R - F_L) * lambda;
        U[idx] = conserVarToRhoPUV(U_, gamma);
    }
}

// swap directions
__global__ void RGRP::kernels::rotate_and_flip_vertically(vec5d *input, vec5d *output, int Nx, int Ny)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x; // column in output
    int y = blockIdx.y * blockDim.y + threadIdx.y; // row in output

    if (x < Nx && y < Ny)
    {
        int outputCol = Ny - 1 - y;
        int outputRow = Nx - 1 - x;
        int outputIdx = outputRow * Ny + outputCol;
        output[outputIdx] = input[y * Nx + x];
        swap(output[outputIdx].data[2], output[outputIdx].data[3]);
    }
}

__global__ void RGRP::kernels::combineResults(vec5d *U1, vec5d *U2, int Nx, int Ny)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index (interface)

    if (i < Nx && j < Ny)
    {
        int idx = j * Nx + i;
        U1[idx] = (U1[idx] + U2[idx]) * 0.5;
    }
}

