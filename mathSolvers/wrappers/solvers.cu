#include "solvers.cuh"
#include "../kernels/ghostCellPresets.cuh"
#include "../kernels/basicFuns.cuh"
using namespace RGRP::kernels;

void RGRP_2D_CUDA::setInitialData(const vector<vector<vec5d>> &U0)
{
    for (auto &row : U0)
    {
        for (auto &u : row)
        {
            _U0_.push_back(u);
        }
    }
    _Nx_ = U0[0].size();
    _Ny_ = U0.size();
    allocateDevice();
}

RGRP_2D_CUDA::~RGRP_2D_CUDA()
{
    if (deviceAllocated)
    {
        freeDevice();
    }
}

void RGRP_2D_CUDA::setRPInitialData(const vec5d &U1, const vec5d &U2, const vec5d &U3, const vec5d &U4, int Nx, int Ny)
{
    _U0_.clear();
    _U0_.resize(Nx * Ny);
    _Nx_ = Nx;
    _Ny_ = Ny;
    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            if (i < Nx / 2 && j < Ny / 2) // 2
            {
                _U0_[j * Nx + i] = U2;
            }
            else if (i < Nx / 2 && j >= Ny / 2) // 3
            {
                _U0_[j * Nx + i] = U3;
            }
            else if (i >= Nx / 2 && j < Ny / 2) // 1
            {
                _U0_[j * Nx + i] = U1;
            }
            else // 4
            {
                _U0_[j * Nx + i] = U4;
            }
        }
    }
    if (deviceAllocated == true)
    {
        freeDevice();
    }
    allocateDevice();
}

void RGRP_2D_CUDA::setSpatialLayout(real width, real height)
{
    _width_ = width;
    _height_ = height;
}

void RGRP_2D_CUDA::setTimeLayout(real endTime, real CFL)
{
    _endTime_ = endTime;
    __CFL__ = CFL;
}

void RGRP_2D_CUDA::setTimeLayout(real endTime)
{
    _endTime_ = endTime;
}

void RGRP_2D_CUDA::setCFL(real CFL)
{
    __CFL__ = CFL;
}

void RGRP_2D_CUDA::setAlpha(real alpha)
{
    _alpha_ = alpha;
}

void RGRP_2D_CUDA::setGamma(real gamma)
{
    _gamma_ = gamma;
}

void RGRP_2D_CUDA::setGhostCellStrategy(ghostCellStrategy strategy)
{
    __strategy__ = strategy;
}

void RGRP_2D_CUDA::allocateDevice()
{
    num_cells = _Nx_ * _Ny_;
    cudaMalloc(&u_dev, num_cells * sizeof(vec5d));
    cudaMalloc(&u_lr_interface_dev, (num_cells + std::max(_Nx_, _Ny_)) * 2 * sizeof(vec5d));
    cudaMalloc(&u_y_dev, num_cells * sizeof(vec5d));
    cudaMalloc(&u_dev_copy, num_cells * sizeof(vec5d));
    cudaMalloc(&u_slope_dev, (num_cells + std::max(_Nx_, _Ny_)) * sizeof(vec5d));
    cudaMalloc(&speed, num_cells * sizeof(real));
    cudaMemset(u_slope_dev, 0, (num_cells + std::max(_Nx_, _Ny_)) * sizeof(vec5d));
    cudaMemcpy(u_dev, _U0_.data(), num_cells * sizeof(vec5d), cudaMemcpyHostToDevice);

    blockSize = dim3(16, 16);
    gridSize = dim3((_Nx_ + blockSize.x - 1) / blockSize.x, (_Ny_ + blockSize.y - 1) / blockSize.y);
    gridSize_y = dim3((_Ny_ + blockSize.y - 1) / blockSize.y, (_Nx_ + blockSize.x - 1) / blockSize.x);
    deviceAllocated = true;
}

void RGRP_2D_CUDA::freeDevice()
{
    cudaFree(u_dev);
    cudaFree(u_lr_interface_dev);
    cudaFree(u_y_dev);
    cudaFree(u_dev_copy);
    cudaFree(u_slope_dev);
    cudaFree(speed);
    deviceAllocated = false;
}

vector<vector<vec5d>> RGRP_2D_CUDA::solve()
{
    int iter = 0;
    int iter_ = 0;
    while (_cTime_ < _endTime_)
    {
        maxPropagtingSpeed<<<gridSize, blockSize>>>(u_dev, speed, _Nx_, _Ny_, _gamma_);
        cudaDeviceSynchronize();
        thrust::device_ptr<real> p_speed = thrust::device_pointer_cast(speed);
        real max_speed = reduce(p_speed, p_speed + num_cells, -std::numeric_limits<real>::infinity(), thrust::maximum<real>());
        real h = std::min(_width_ / real(_Nx_), _height_ / real(_Ny_));
        real dt = std::min(__CFL__ * h / max_speed, _endTime_ - _cTime_);
        if (dt < 1e-20)
        {
            break;
        }

        cudaMemcpy(u_dev_copy, u_dev, num_cells * sizeof(vec5d), cudaMemcpyDeviceToDevice);

        MUSCL<<<gridSize, blockSize>>>(u_dev, u_slope_dev, _Nx_, _Ny_, _alpha_, _width_, u_lr_interface_dev);
        cudaDeviceSynchronize();
        rel_RP_posi_fix<<<gridSize, blockSize>>>(u_dev, u_slope_dev, u_lr_interface_dev, _Nx_, _Ny_, _gamma_);
        cudaDeviceSynchronize();
        forward_x_dir<<<gridSize, blockSize>>>(u_dev, u_lr_interface_dev, _Nx_, _Ny_, _width_, dt, _gamma_);
        cudaDeviceSynchronize();
        rotate_and_flip_vertically<<<gridSize, blockSize>>>(u_dev, u_y_dev, _Nx_, _Ny_);
        cudaDeviceSynchronize();
        MUSCL<<<gridSize_y, blockSize>>>(u_y_dev, u_slope_dev, _Ny_, _Nx_, _alpha_, _height_, u_lr_interface_dev);
        cudaDeviceSynchronize();
        rel_RP_posi_fix<<<gridSize_y, blockSize>>>(u_y_dev, u_slope_dev, u_lr_interface_dev, _Ny_, _Nx_, _gamma_);
        cudaDeviceSynchronize();
        forward_x_dir<<<gridSize_y, blockSize>>>(u_y_dev, u_lr_interface_dev, _Ny_, _Nx_, _height_, dt, _gamma_);
        cudaDeviceSynchronize();
        rotate_and_flip_vertically<<<gridSize_y, blockSize>>>(u_y_dev, u_dev, _Ny_, _Nx_);
        cudaDeviceSynchronize();

        rotate_and_flip_vertically<<<gridSize, blockSize>>>(u_dev_copy, u_y_dev, _Nx_, _Ny_);
        cudaDeviceSynchronize();
        MUSCL<<<gridSize_y, blockSize>>>(u_y_dev, u_slope_dev, _Ny_, _Nx_, _alpha_, _height_, u_lr_interface_dev);
        cudaDeviceSynchronize();
        rel_RP_posi_fix<<<gridSize_y, blockSize>>>(u_y_dev, u_slope_dev, u_lr_interface_dev, _Ny_, _Nx_, _gamma_);
        cudaDeviceSynchronize();
        forward_x_dir<<<gridSize_y, blockSize>>>(u_y_dev, u_lr_interface_dev, _Ny_, _Nx_, _height_, dt, _gamma_);
        cudaDeviceSynchronize();
        rotate_and_flip_vertically<<<gridSize_y, blockSize>>>(u_y_dev, u_dev_copy, _Ny_, _Nx_);
        cudaDeviceSynchronize();
        MUSCL<<<gridSize, blockSize>>>(u_dev_copy, u_slope_dev, _Nx_, _Ny_, _alpha_, _width_, u_lr_interface_dev);
        cudaDeviceSynchronize();
        rel_RP_posi_fix<<<gridSize, blockSize>>>(u_dev_copy, u_slope_dev, u_lr_interface_dev, _Nx_, _Ny_, _gamma_);
        cudaDeviceSynchronize();
        forward_x_dir<<<gridSize, blockSize>>>(u_dev_copy, u_lr_interface_dev, _Nx_, _Ny_, _width_, dt, _gamma_);
        cudaDeviceSynchronize();

        combineResults<<<gridSize, blockSize>>>(u_dev, u_dev_copy, _Nx_, _Ny_);

        _cTime_ += dt;
        iter++;
        iter_++;

        if (__strategy__ == ghostCellStrategy::outflow)
        {
            set_ghost_cells<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
        }
        else if (__strategy__ == ghostCellStrategy::jet)
        {
            set_ghost_cells_jet<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
        }
        else if (__strategy__ == ghostCellStrategy::jet800)
        {
            set_ghost_cells_jet800<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
        }
        else if (__strategy__ == ghostCellStrategy::RMInstability)
        {
            set_ghost_cells_RMInstability<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
        }
        else
        {
            set_ghost_cells_2Mach<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
        }
        cudaDeviceSynchronize();
    }

    cudaMemcpy(_U0_.data(), u_dev, num_cells * sizeof(vec5d), cudaMemcpyDeviceToHost);
    _U_ = flattened2vectorized(_U0_);
    return _U_;
}

void RGRP_2D_CUDA::iterateOnce()
{
    maxPropagtingSpeed<<<gridSize, blockSize>>>(u_dev, speed, _Nx_, _Ny_, _gamma_);
    cudaDeviceSynchronize();
    thrust::device_ptr<real> p_speed = thrust::device_pointer_cast(speed);
    real max_speed = reduce(p_speed, p_speed + num_cells, -std::numeric_limits<real>::infinity(), thrust::maximum<real>());
    real h = std::min(_width_ / real(_Nx_), _height_ / real(_Ny_));
    real dt = std::min(__CFL__ * h / max_speed, _endTime_ - _cTime_);
    if (dt < 1e-20)
    {
        notCompleted = false;
        return;
    }

    cudaMemcpy(u_dev_copy, u_dev, num_cells * sizeof(vec5d), cudaMemcpyDeviceToDevice);

    MUSCL<<<gridSize, blockSize>>>(u_dev, u_slope_dev, _Nx_, _Ny_, _alpha_, _width_, u_lr_interface_dev);
    cudaDeviceSynchronize();
    rel_RP_posi_fix<<<gridSize, blockSize>>>(u_dev, u_slope_dev, u_lr_interface_dev, _Nx_, _Ny_, _gamma_);
    cudaDeviceSynchronize();
    forward_x_dir<<<gridSize, blockSize>>>(u_dev, u_lr_interface_dev, _Nx_, _Ny_, _width_, dt, _gamma_);
    cudaDeviceSynchronize();
    rotate_and_flip_vertically<<<gridSize, blockSize>>>(u_dev, u_y_dev, _Nx_, _Ny_);
    cudaDeviceSynchronize();
    MUSCL<<<gridSize_y, blockSize>>>(u_y_dev, u_slope_dev, _Ny_, _Nx_, _alpha_, _height_, u_lr_interface_dev);
    cudaDeviceSynchronize();
    rel_RP_posi_fix<<<gridSize_y, blockSize>>>(u_y_dev, u_slope_dev, u_lr_interface_dev, _Ny_, _Nx_, _gamma_);
    cudaDeviceSynchronize();
    forward_x_dir<<<gridSize_y, blockSize>>>(u_y_dev, u_lr_interface_dev, _Ny_, _Nx_, _height_, dt, _gamma_);
    cudaDeviceSynchronize();
    rotate_and_flip_vertically<<<gridSize_y, blockSize>>>(u_y_dev, u_dev, _Ny_, _Nx_);
    cudaDeviceSynchronize();

    rotate_and_flip_vertically<<<gridSize, blockSize>>>(u_dev_copy, u_y_dev, _Nx_, _Ny_);
    cudaDeviceSynchronize();
    MUSCL<<<gridSize_y, blockSize>>>(u_y_dev, u_slope_dev, _Ny_, _Nx_, _alpha_, _height_, u_lr_interface_dev);
    cudaDeviceSynchronize();
    rel_RP_posi_fix<<<gridSize_y, blockSize>>>(u_y_dev, u_slope_dev, u_lr_interface_dev, _Ny_, _Nx_, _gamma_);
    cudaDeviceSynchronize();
    forward_x_dir<<<gridSize_y, blockSize>>>(u_y_dev, u_lr_interface_dev, _Ny_, _Nx_, _height_, dt, _gamma_);
    cudaDeviceSynchronize();
    rotate_and_flip_vertically<<<gridSize_y, blockSize>>>(u_y_dev, u_dev_copy, _Ny_, _Nx_);
    cudaDeviceSynchronize();
    MUSCL<<<gridSize, blockSize>>>(u_dev_copy, u_slope_dev, _Nx_, _Ny_, _alpha_, _width_, u_lr_interface_dev);
    cudaDeviceSynchronize();
    rel_RP_posi_fix<<<gridSize, blockSize>>>(u_dev_copy, u_slope_dev, u_lr_interface_dev, _Nx_, _Ny_, _gamma_);
    cudaDeviceSynchronize();
    forward_x_dir<<<gridSize, blockSize>>>(u_dev_copy, u_lr_interface_dev, _Nx_, _Ny_, _width_, dt, _gamma_);
    cudaDeviceSynchronize();

    combineResults<<<gridSize, blockSize>>>(u_dev, u_dev_copy, _Nx_, _Ny_);

    _cTime_ += dt;
    iter++;

    if (__strategy__ == ghostCellStrategy::outflow)
    {
        set_ghost_cells<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
    }
    else if (__strategy__ == ghostCellStrategy::jet)
    {
        set_ghost_cells_jet<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
    }
    else if (__strategy__ == ghostCellStrategy::jet800)
    {
        set_ghost_cells_jet800<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
    }
    else if (__strategy__ == ghostCellStrategy::RMInstability)
    {
        set_ghost_cells_RMInstability<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
    }
    else
    {
        set_ghost_cells_2Mach<<<gridSize, blockSize>>>(u_dev, _Nx_, _Ny_, _width_, _height_, _cTime_);
    }
    cudaDeviceSynchronize();
}

void RGRP_2D_CUDA::copyFromDevice()
{
    cudaMemcpy(_U0_.data(), u_dev, num_cells * sizeof(vec5d), cudaMemcpyDeviceToHost);
    _U_ = flattened2vectorized(_U0_);
}

vector<vector<vec5d>> RGRP_2D_CUDA::flattened2vectorized(const vector<vec5d> &flattened)
{
    vector<vector<vec5d>> ret;
    for (int j = 0; j < _Ny_; ++j)
    {
        vector<vec5d> temp;
        for (int i = 0; i < _Nx_; ++i)
        {
            temp.push_back(flattened[j * _Nx_ + i]);
        }
        ret.push_back(temp);
    }
    return ret;
}