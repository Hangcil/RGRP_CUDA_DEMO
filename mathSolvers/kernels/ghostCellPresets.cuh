#pragma once
#include "vec.cuh"
using RGRP::kernels::vec5d;

namespace RGRP::kernels
{
    // ghost cell setup for the real Mach reflection problem
    __global__ void set_ghost_cells_2Mach(vec5d *U, int Nx, int Ny, real width, real height, real t);

    // ghost cell setup for the jet problem
    __global__ void set_ghost_cells_jet(vec5d *U, int Nx, int Ny, real width, real height, real t);

    // ghost cell setup for the jet 800 problem
    __global__ void set_ghost_cells_jet800(vec5d *U, int Nx, int Ny, real width, real height, real t);

    __global__ void set_ghost_cells(vec5d *U, int Nx, int Ny, real width, real height, real t);

    // ghost cell setup for the Richtmyer-Meshkov instability problem
    __global__ void set_ghost_cells_RMInstability(vec5d *U, int Nx, int Ny, real width, real height, real t);
}