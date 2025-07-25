#include "ghostCellPresets.cuh"

// ghost cell setup for the real Mach reflection problem
__global__ void RGRP::kernels::set_ghost_cells_2Mach(vec5d *U, int Nx, int Ny, real width, real height, real t)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index (interface)

    int idx = j * Nx + i;
    if (i < Nx && j < Ny)
    {
        if (j == 0 || j == 1)
        {
            real temp = real(Ny) / sqrt(3.0);
            real dist = real(Ny) * 20 * t / sqrt(3.0);
            if (i <= Nx / 24 + int(temp) + int(dist))
            {
                U[idx] = {8, 116.5, 4.125 * sqrt(3.0), -4.125, 0};
            }
            else
            {
                U[idx] = {1.4, 1, 0, 0, 0};
            }
        }
        else if (j == Ny - 1)
        {
            if (i <= Nx / 24)
            {
                U[idx] = U[idx - 2 * Nx];
            }
            else
            {
                U[idx] = U[idx - 3 * Nx];
                U[idx](3) = -U[idx](3);
            }
        }
        else if (j == Ny - 2)
        {
            if (i <= Nx / 24)
            {
                U[idx] = U[idx - Nx];
            }
            else
            {
                U[idx] = U[idx - Nx];
                U[idx](3) = -U[idx](3);
            }
        }

        if (i == 0 || i == 1)
        {
            U[idx] = {8, 116.5, 4.125 * sqrt(3.0), -4.125, 0};
        }
        else if (i == Nx - 2 || i == Nx - 1)
        {
            U[idx] = {1.4, 1, 0, 0, 0};
        }
    }
}

// ghost cell setup for outflow boundary conditions
__global__ void RGRP::kernels::set_ghost_cells(vec5d *U, int Nx, int Ny, real width, real height, real t)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index (interface)

    int idx = j * Nx + i;
    if (j == 0)
    {
        U[idx] = U[idx + 2 * Nx];
    }
    else if (j == 1)
    {
        U[idx] = U[idx + Nx];
    }
    else if (j == Ny - 2)
    {
        U[idx] = U[idx - Nx];
    }
    else if (j == Ny - 1)
    {
        U[idx] = U[idx - 2 * Nx];
    }

    if (i == 0)
    {
        U[idx] = U[idx + 2];
    }
    else if (i == 1)
    {
        U[idx] = U[idx + 1];
    }
    else if (i == Nx - 2)
    {
        U[idx] = U[idx - 1];
    }
    else if (i == Nx - 1)
    {
        U[idx] = U[idx - 2];
    }
}

// ghost cell setup for the jet problem
__global__ void RGRP::kernels::set_ghost_cells_jet(vec5d *U, int Nx, int Ny, real width, real height, real t)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index (interface)

    if (i < Nx && j < Ny)
    {
        int idx = j * Nx + i;
        if (j == 0)
        {
            U[idx] = U[idx + 2 * Nx];
        }
        else if (j == 1)
        {
            U[idx] = U[idx + Nx];
        }
        else if (j == Ny - 2)
        {
            U[idx] = U[idx - Nx];
        }
        else if (j == Ny - 1)
        {
            U[idx] = U[idx - 2 * Nx];
        }

        if (i == 0 || i == 1)
        {
            if (j >= 2 * Ny / 5 + 1 && j <= 3 * Ny / 5)
            {
                U[idx] = vec5d(5.0, 0.4127, 800.0, 0.0, 0.0);
            }
            else
            {
                U[idx] = {0.5, 0.4127, 0.0, 0.0, 0.0};
            }
        }
        else if (i == Nx - 2)
        {
            U[idx] = U[idx - 1];
        }
        else if (i == Nx - 1)
        {
            U[idx] = U[idx - 2];
        }
    }
}

// ghost cell setup for the jet 800 problem
__global__ void RGRP::kernels::set_ghost_cells_jet800(vec5d *U, int Nx, int Ny, real width, real height, real t)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index (interface)

    if (i < Nx && j < Ny)
    {
        int idx = j * Nx + i;
        if (j == 0)
        {
            U[idx] = U[idx + 2 * Nx];
        }
        else if (j == 1)
        {
            U[idx] = U[idx + Nx];
        }
        else if (j == Ny - 2)
        {
            U[idx] = U[idx - Nx];
            U[idx](3) = -U[idx](3);
        }
        else if (j == Ny - 1)
        {
            U[idx] = U[idx - 3 * Nx];
            U[idx](3) = -U[idx](3);
        }

        if (i == 0 || i == 1)
        {
            if (j >= 9 * Ny / 10)
            {
                U[idx] = vec5d(1.4, 1, 800, 0, 0);
            }
            else
            {
                U[idx] = {0.14, 1, 0, 0, 0};
            }
        }
        else if (i == Nx - 2)
        {
            U[idx] = U[idx - 1];
        }
        else if (i == Nx - 1)
        {
            U[idx] = U[idx - 2];
        }
    }
}

// ghost cell setup for the Richtmyer-Meshkov instability problem
__global__ void RGRP::kernels::set_ghost_cells_RMInstability(vec5d *U, int Nx, int Ny, real width, real height, real t)
{
    int j = blockIdx.y * blockDim.y + threadIdx.y; // y-index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // x-index (interface)

    if (i < Nx && j < Ny)
    {
        int idx = j * Nx + i;
        if (j == 0 || j == 1)
        {
            U[idx] = {8.0 / 3.0, 4.5, 0, 0, 0};
        }
        else if (j == Ny - 1 || j == Ny - 2)
        {
            U[idx] = {0.25, 1, 0, 0, 0};
        }

        if (i == 0)
        {
            U[idx] = U[idx + 3];
            U[idx](2) = -U[idx](2);
        }
        else if (i == 1)
        {
            U[idx] = U[idx + 1];
            U[idx](2) = -U[idx](2);
        }
        else if (i == Nx - 2)
        {
            U[idx] = U[idx - 1];
            U[idx](2) = -U[idx](2);
        }
        else if (i == Nx - 1)
        {
            U[idx] = U[idx - 3];
            U[idx](2) = -U[idx](2);
        }
    }
}
