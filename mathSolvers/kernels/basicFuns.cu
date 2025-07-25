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

        real rho_l = U_lr_interface[2 * idx](0);
        real rho_r = U_lr_interface[2 * idx + 1](0);
        real pi_l = U_lr_interface[2 * idx](1);
        real pi_r = U_lr_interface[2 * idx + 1](1);
        real u_l = U_lr_interface[2 * idx](2);
        real u_r = U_lr_interface[2 * idx + 1](2);
        real v_l = U_lr_interface[2 * idx](3);
        real v_r = U_lr_interface[2 * idx + 1](3);
        real e_l = pi_l / rho_l / (gamma - 1.0);
        real e_r = pi_r / rho_r / (gamma - 1.0);
        real c_l = 0.0;
        real c_r = 0.0;
        real rho_lx = Us[idxL](0);
        real rho_rx = Us[idx](0);
        real pi_lx = Us[idxL](1);
        real pi_rx = Us[idx](1);
        real u_lx = Us[idxL](2);
        real u_rx = Us[idx](2);
        real v_lx = Us[idxL](3);
        real v_rx = Us[idx](3);
        real e_lx = pi_lx / rho_l / (gamma - 1.0) - pi_l * rho_lx / rho_l / rho_l / (gamma - 1.0);
        real e_rx = pi_rx / rho_r / (gamma - 1.0) - pi_r * rho_rx / rho_r / rho_r / (gamma - 1.0);
        real c_lx = 0.5 * sqrt(gamma / pi_l / rho_l) * (pi_l * rho_lx + rho_l * pi_lx);
        real c_rx = 0.5 * sqrt(gamma / pi_r / rho_r) * (pi_r * rho_rx + rho_r * pi_rx);
        real rho_t = 0.0, pi_t = 0.0, u_t = 0.0, v_t = 0.0, e_t = 0.0;

        real a_l = sqrt(gamma * pi_l / rho_l), a_r = sqrt(gamma * pi_r / rho_r);
        real kappa = (gamma + 1.0) / 2.0;
        if (pi_r >= pi_l)
        {
            /*real c_l_*/ a_l = rho_l * (a_l + kappa * max(0.0, (pi_r - pi_l) / rho_r / a_r + u_l - u_r));
            c_l = /*c_l_*/ a_l;
            c_r = rho_r * (a_r + kappa * max(0.0, (pi_l - pi_r) / /* c_l_ */ a_l + u_l - u_r));
        }
        else
        {
            /* real c_r_ */ a_r = rho_r * (a_r + kappa * max(0.0, (pi_l - pi_r) / rho_l / a_l + u_l - u_r));
            /* real c_l_ */ a_l = rho_l * (a_l + kappa * max(0.0, (pi_r - pi_l) / /* c_r_ */ a_r + u_l - u_r));
            c_l = /* c_l_ */ a_l;
            c_r = /* c_r_ */ a_r;
        }

        if (u_l - c_l / rho_l >= 0.0)
        {
            rho_t = -rho_l * u_lx - u_l * rho_lx;
            u_t = -u_l * u_lx - pi_lx / rho_l;
            e_t = -u_l * e_lx - pi_l / rho_l * u_lx;
            v_t = -u_l * v_lx;
            pi_t = -u_l * pi_lx - c_l * c_l / rho_l * u_lx;
            U_lr_interface[2 * idx] = vec5d{rho_l, pi_l, u_l, v_l, e_l};
        }
        else if (u_r + c_r / rho_r <= 0.0)
        {
            rho_t = -rho_r * u_rx - u_r * rho_rx;
            u_t = -u_r * u_rx - pi_rx / rho_r;
            e_t = -u_r * e_rx - pi_r / rho_r * u_rx;
            v_t = -u_r * v_rx;
            pi_t = -u_r * pi_rx - c_r * c_r / rho_r * u_rx;
            U_lr_interface[2 * idx] = vec5d{rho_r, pi_r, u_r, v_r, e_r};
        }
        else
        {
            real u_m = (u_l * c_l + u_r * c_r + pi_l - pi_r) / (c_l + c_r);
            real pi_m = ((u_l - u_r) * c_l * c_r + c_l * pi_r + c_r * pi_l) / (c_l + c_r);
            real k_l = -pi_lx - c_l * u_lx + 0.5 * (u_m - u_l) * c_lx;
            real k_r = -pi_rx + c_r * u_rx - 0.5 * (u_m - u_r) * c_rx;
            real d_l = -2.0 * c_l * c_l / rho_l * u_lx - 2.0 * c_l / rho_l * pi_lx + c_l / rho_l * (u_m - u_l) * c_lx;
            real d_r = 2.0 * c_r * c_r / rho_r * u_rx - 2.0 * c_r / rho_r * pi_rx + c_r / rho_r * (u_r - u_m) * c_rx;
            real Dpi_Dt = (c_r * d_l - c_l * d_r) / 2.0 / (c_l + c_r);
            if (u_m <= 0.0)
            {
                real rho_mr = c_r / (u_r - u_m + c_r / rho_r);
                real e_mr = e_r + (pi_m - pi_r) * (pi_m + pi_r) / 2.0 / c_r / c_r;
                rho_t = -u_m * rho_mr * rho_mr * rho_mr / (rho_r * c_r * c_r) * (c_r * c_r / rho_r / rho_r * rho_rx - c_r * u_rx + 1.5 * (u_m - u_r) * c_rx) + rho_mr * rho_mr * rho_mr / (c_r * c_r * c_r) * (u_m + c_r / rho_mr) * Dpi_Dt;
                pi_t = 1.0 / (c_l + c_r) * (c_l * c_r / rho_l * (1 + u_m * rho_mr / c_r) * k_l - c_l * c_r / rho_r * (1 - u_m * rho_mr / c_l) * k_r);
                u_t = 1.0 / (c_l + c_r) * (c_l / rho_l * (1 + u_m * rho_mr / c_r) * k_l + c_r / rho_r * (1 - u_m * rho_mr * c_l / c_r / c_r) * k_r);
                e_t = -u_m * rho_mr / rho_r / c_r / c_r / c_r * (-pi_r * c_r * pi_rx + c_r * c_r * c_r * e_rx + 2.0 * (e_r - e_mr) * c_r * c_r * c_rx) + pi_m / c_r / c_r * pi_t;
                v_t = -u_m * rho_mr / rho_r * v_rx;
                U_lr_interface[2 * idx] = vec5d{rho_mr, pi_m, u_m, v_r, e_mr};
            }
            else
            {
                real rho_ml = c_l / (u_m - u_l + c_l / rho_l);
                real e_ml = e_l + (pi_m - pi_l) * (pi_m + pi_l) / 2.0 / c_l / c_l;
                rho_t = -u_m * rho_ml * rho_ml * rho_ml / (rho_l * c_l * c_l) * (c_l * c_l / rho_l / rho_l * rho_lx + c_l * u_lx + 1.5 * (u_m - u_l) * c_lx) - rho_ml * rho_ml * rho_ml / (c_l * c_l * c_l) * (u_m - c_l / rho_ml) * Dpi_Dt;
                pi_t = 1.0 / (c_l + c_r) * (c_l * c_r / rho_l * (1 + u_m * rho_ml / c_r) * k_l - c_l * c_r / rho_r * (1 - u_m * rho_ml / c_l) * k_r);
                u_t = 1.0 / (c_l + c_r) * (c_l / rho_l * (1 + u_m * rho_ml * c_r / c_l / c_l) * k_l + c_r / rho_r * (1 - u_m * rho_ml / c_l) * k_r);
                e_t = -u_m * rho_ml / rho_l / c_l / c_l / c_l * (-pi_l * c_l * pi_lx + c_l * c_l * c_l * e_lx + 2.0 * (e_l - e_ml) * c_l * c_l * c_lx) + pi_m / c_l / c_l * pi_t;
                v_t = -u_m * rho_ml / rho_l * v_lx;
                U_lr_interface[2 * idx] = vec5d{rho_ml, pi_m, u_m, v_l, e_ml};
            }
        }
        U_lr_interface[2 * idx + 1] = vec5d(rho_t, pi_t, u_t, v_t, e_t);
    }
}

__global__ void RGRP::kernels::GRP_l(vec5d *U, vec5d *Us, vec5d *U_lr_interface, int Nx, int Ny, real gamma)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index

    if (i < Nx - 1 && j < Ny && i > 1)
    {
        int idxL = j * Nx + i - 1; // left cell
        int idx = j * Nx + i;      // cell i(left interface)

        real rho_l = U_lr_interface[2 * idx](0);
        real rho_r = U_lr_interface[2 * idx + 1](0);
        real pi_l = U_lr_interface[2 * idx](1);
        real pi_r = U_lr_interface[2 * idx + 1](1);
        real u_l = U_lr_interface[2 * idx](2);
        real u_r = U_lr_interface[2 * idx + 1](2);
        real v_l = U_lr_interface[2 * idx](3);
        real v_r = U_lr_interface[2 * idx + 1](3);
        real e_l = pi_l / rho_l / (gamma - 1.0);
        real e_r = pi_r / rho_r / (gamma - 1.0);
        real c_l = 0.0;
        real c_r = 0.0;
        real rho_lx = Us[idxL](0);
        real rho_rx = Us[idx](0);
        real pi_lx = Us[idxL](1);
        real pi_rx = Us[idx](1);
        real u_lx = Us[idxL](2);
        real u_rx = Us[idx](2);
        real v_lx = Us[idxL](3);
        real v_rx = Us[idx](3);
        real e_lx = pi_lx / rho_l / (gamma - 1.0) - pi_l * rho_lx / rho_l / rho_l / (gamma - 1.0);
        real e_rx = pi_rx / rho_r / (gamma - 1.0) - pi_r * rho_rx / rho_r / rho_r / (gamma - 1.0);
        real c_lx = 0.5 * sqrt(gamma / pi_l / rho_l) * (pi_l * rho_lx + rho_l * pi_lx);
        real c_rx = 0.5 * sqrt(gamma / pi_r / rho_r) * (pi_r * rho_rx + rho_r * pi_rx);
        real rho_t = 0.0, pi_t = 0.0, u_t = 0.0, v_t = 0.0, e_t = 0.0;

        real a_l = sqrt(gamma * pi_l / rho_l), a_r = sqrt(gamma * pi_r / rho_r);
        real kappa = (gamma + 1.0) / 2.0;
        if (pi_r >= pi_l)
        {
            /*real c_l_*/ a_l = rho_l * (a_l + kappa * max(0.0, (pi_r - pi_l) / rho_r / a_r + u_l - u_r));
            c_l = /*c_l_*/ a_l;
            c_r = rho_r * (a_r + kappa * max(0.0, (pi_l - pi_r) / /* c_l_ */ a_l + u_l - u_r));
        }
        else
        {
            /* real c_r_ */ a_r = rho_r * (a_r + kappa * max(0.0, (pi_l - pi_r) / rho_l / a_l + u_l - u_r));
            /* real c_l_ */ a_l = rho_l * (a_l + kappa * max(0.0, (pi_r - pi_l) / /* c_r_ */ a_r + u_l - u_r));
            c_l = /* c_l_ */ a_l;
            c_r = /* c_r_ */ a_r;
        }

        if (u_l - c_l / rho_l >= 0.0)
        {
            rho_t = -rho_l * u_lx - u_l * rho_lx;
            u_t = -u_l * u_lx - pi_lx / rho_l;
            e_t = -u_l * e_lx - pi_l / rho_l * u_lx;
            v_t = -u_l * v_lx;
            pi_t = -u_l * pi_lx - c_l * c_l / rho_l * u_lx;
            U_lr_interface[2 * idx] = vec5d{rho_l, pi_l, u_l, v_l, e_l};
        }
        else if (u_r + c_r / rho_r <= 0.0)
        {
            rho_t = -rho_r * u_rx - u_r * rho_rx;
            u_t = -u_r * u_rx - pi_rx / rho_r;
            e_t = -u_r * e_rx - pi_r / rho_r * u_rx;
            v_t = -u_r * v_rx;
            pi_t = -u_r * pi_rx - c_r * c_r / rho_r * u_rx;
            U_lr_interface[2 * idx] = vec5d{rho_r, pi_r, u_r, v_r, e_r};
        }
        else
        {
            real u_m = (u_l * c_l + u_r * c_r + pi_l - pi_r) / (c_l + c_r);
            real pi_m = ((u_l - u_r) * c_l * c_r + c_l * pi_r + c_r * pi_l) / (c_l + c_r);
            real k_l = -pi_lx - c_l * u_lx + 0.5 * (u_m - u_l) * c_lx;
            real k_r = -pi_rx + c_r * u_rx - 0.5 * (u_m - u_r) * c_rx;
            real d_l = -2.0 * c_l * c_l / rho_l * u_lx - 2.0 * c_l / rho_l * pi_lx + c_l / rho_l * (u_m - u_l) * c_lx;
            real d_r = 2.0 * c_r * c_r / rho_r * u_rx - 2.0 * c_r / rho_r * pi_rx + c_r / rho_r * (u_r - u_m) * c_rx;
            real Dpi_Dt = (c_r * d_l - c_l * d_r) / 2.0 / (c_l + c_r);
            if (u_m <= 0.0)
            {
                real rho_mr = c_r / (u_r - u_m + c_r / rho_r);
                real e_mr = e_r + (pi_m - pi_r) * (pi_m + pi_r) / 2.0 / c_r / c_r;
                rho_t = -u_m * rho_mr * rho_mr * rho_mr / (rho_r * c_r * c_r) * (c_r * c_r / rho_r / rho_r * rho_rx - c_r * u_rx + 1.5 * (u_m - u_r) * c_rx) + rho_mr * rho_mr * rho_mr / (c_r * c_r * c_r) * (u_m + c_r / rho_mr) * Dpi_Dt;
                pi_t = 1.0 / (c_l + c_r) * (c_l * c_r / rho_l * (1 + u_m * rho_mr / c_r) * k_l - c_l * c_r / rho_r * (1 - u_m * rho_mr / c_l) * k_r);
                u_t = 1.0 / (c_l + c_r) * (c_l / rho_l * (1 + u_m * rho_mr / c_r) * k_l + c_r / rho_r * (1 - u_m * rho_mr * c_l / c_r / c_r) * k_r);
                e_t = -u_m * rho_mr / rho_r / c_r / c_r / c_r * (-pi_r * c_r * pi_rx + c_r * c_r * c_r * e_rx + 2.0 * (e_r - e_mr) * c_r * c_r * c_rx) + pi_m / c_r / c_r * pi_t;
                v_t = -u_m * rho_mr / rho_r * v_rx;
                U_lr_interface[2 * idx] = vec5d{rho_mr, pi_m, u_m, v_r, e_mr};
            }
            else
            {
                real rho_ml = c_l / (u_m - u_l + c_l / rho_l);
                real e_ml = e_l + (pi_m - pi_l) * (pi_m + pi_l) / 2.0 / c_l / c_l;
                rho_t = -u_m * rho_ml * rho_ml * rho_ml / (rho_l * c_l * c_l) * (c_l * c_l / rho_l / rho_l * rho_lx + c_l * u_lx + 1.5 * (u_m - u_l) * c_lx) - rho_ml * rho_ml * rho_ml / (c_l * c_l * c_l) * (u_m - c_l / rho_ml) * Dpi_Dt;
                pi_t = 1.0 / (c_l + c_r) * (c_l * c_r / rho_l * (1 + u_m * rho_ml / c_r) * k_l - c_l * c_r / rho_r * (1 - u_m * rho_ml / c_l) * k_r);
                u_t = 1.0 / (c_l + c_r) * (c_l / rho_l * (1 + u_m * rho_ml * c_r / c_l / c_l) * k_l + c_r / rho_r * (1 - u_m * rho_ml / c_l) * k_r);
                e_t = -u_m * rho_ml / rho_l / c_l / c_l / c_l * (-pi_l * c_l * pi_lx + c_l * c_l * c_l * e_lx + 2.0 * (e_l - e_ml) * c_l * c_l * c_lx) + pi_m / c_l / c_l * pi_t;
                v_t = -u_m * rho_ml / rho_l * v_lx;
                U_lr_interface[2 * idx] = vec5d{rho_ml, pi_m, u_m, v_l, e_ml};
            }
        }
        U_lr_interface[2 * idx + 1] = vec5d(rho_t, pi_t, u_t, v_t, e_t);
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
