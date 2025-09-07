#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "../kernels/vec.cuh"
using real = double;

using std::ofstream;
using std::string;
using std::vector;
using RGRP::kernels::vec5d;
enum class ghostCellStrategy
{
    outflow,
    reflective,
    jet,
    jet800,
    doubleMach,
    RMInstability
};

class RGRP_2D_CUDA
{
public:
    void setInitialData(const vector<vector<vec5d>> &U0);
    RGRP_2D_CUDA() {}
    ~RGRP_2D_CUDA();
    void setRPInitialData(const vec5d &U1, const vec5d &U2, const vec5d &U3, const vec5d &U4, int Nx, int Ny);
    void setSpatialLayout(real width, real height);
    void setTimeLayout(real endTime, real CFL);
    void setTimeLayout(real endTime);
    void setCFL(real CFL);
    void setAlpha(real alpha);
    void setGamma(real gamma);
    void setGhostCellStrategy(ghostCellStrategy strategy);
    void allocateDevice();
    void freeDevice();
    void iterateOnce();
    void copyFromDevice();


    int iter = 0;
    bool notCompleted = true;
    real _cTime_ = 0.0;
    vector<vector<vec5d>> _U_;

protected:
    vector<vec5d> _U0_;
    real _width_ = 1.0, _height_ = 1.0, _endTime_ = 0.0, __CFL__ = 0.5;
    real _alpha_ = 1.9, _gamma_ = 1.4;
    int _Nx_ = 100, _Ny_ = 100, num_cells = 0;
    ghostCellStrategy __strategy__ = ghostCellStrategy::outflow;

    vec5d *u_dev;
    vec5d *u_dev_copy;
    vec5d *u_lr_interface_dev;
    vec5d *u_y_dev;
    vec5d *u_slope_dev;
    real *speed;
    dim3 blockSize;
    dim3 gridSize;
    dim3 gridSize_y;
    bool deviceAllocated = false;

    vector<vector<vec5d>> flattened2vectorized(const vector<vec5d> &flattened);
};

