#include "kernelLauncher.h"

kernelLauncher::kernelLauncher(QObject *parent)
{
}

QImage kernelLauncher::convertToQImage(const vector<vector<vec5d>> &data)
{
    const size_t rows = data.size();
    const size_t cols = data[0].size();
    double min_val = data[0][0](0);
    double max_val = data[0][0](0);

    for (const auto &row : data)
    {
        for (auto &val : row)
        {
            min_val = std::min(min_val, val(0));
            max_val = std::max(max_val, val(0));
        }
    }
    if (max_val == min_val)
    {
        min_val = max_val - 1.0;
    }

    QImage heatmap(cols, rows, QImage::Format_RGB32);
    const int colormapSize = RdBuColormap.size();
    for (size_t y = 0; y < rows; ++y)
    {
        for (size_t x = 0; x < cols; ++x)
        {
            double normalized = (data[y][x](0) - min_val) / (max_val - min_val);
            normalized = std::clamp(normalized, 0.0, 1.0);
            int index = static_cast<int>(normalized * (colormapSize - 1));
            heatmap.setPixel(x, y, RdBuColormap[index]);
        }
    }
    return heatmap;
}

void kernelLauncher::writeToTxt(const QString &dir, const vector<vector<vec5d>> &data)
{
    QStringList dirs(4);
    dirs[0] = dir + "/rho.txt";
    dirs[1] = dir + "/p.txt";
    dirs[2] = dir + "/u.txt";
    dirs[3] = dir + "/v.txt";
#pragma omp parallel for
    for (int i = 0; i < 4; i++)
    {
        ofstream of(dirs[i].toStdString());
        if (of.is_open())
        {
            for (const auto &row : data)
            {
                for (const auto &elem : row)
                {
                    of << elem[i] << " ";
                }
                of << "\n";
            }
        }
        of.close();
    }
}

void kernelLauncher::reset()
{
    T = 0.0;
    emit progressChanged();
}

Q_INVOKABLE void kernelLauncher::launchWith(double CFL, double endT, double gamma, double alpha, int Nx, int Ny, int testID)
{
    vector<vector<vec5d>> U0;
    ghostCellStrategy strategy;
    double width = 1.0, height = 1.0;
    T = endT;
    switch (testID)
    {
    case 1:
        U0 = initialize_4shocks(Nx);
        strategy = ghostCellStrategy::outflow;
        break;
    case 2:
        U0 = initialize_4cds(Nx);
        strategy = ghostCellStrategy::outflow;
        break;
    case 3:
        U0 = initialize_jet800(Ny);
        strategy = ghostCellStrategy::jet800;
        width = 1.5;
        height = 0.5;
        break;
    case 4:
        U0 = initialize_2Mach(Ny);
        strategy = ghostCellStrategy::doubleMach;
        width = 4.0;
        height = 1.0;
        break;
    default:
        break;
    }
    RGRP_2D_CUDA solver;
    solver.setInitialData(U0);
    solver.setSpatialLayout(width, height);
    solver.setTimeLayout(endT, CFL);
    solver.setGhostCellStrategy(strategy);
    solver.setAlpha(alpha);
    solver.setGamma(gamma);
    double period = T / double(incrementCount);
    double nextCheckpoint = period;
    while (solver._cTime_ < T)
    {
        solver.iterateOnce();
        if (solver._cTime_ >= nextCheckpoint)
        {
            emit progressChanged();
            nextCheckpoint += period;
        }
    }
    emit simulationFinished();

    solver.copyFromDevice();
    plot = convertToQImage(solver._U_);
    plot.save(outputDir + "/data.png", "PNG");
    emit plotFinished(outputDir + "/data.png");
}

std::vector<std::vector<vec5d>> kernelLauncher::initialize_4shocks(int N)
{
    vec5d U1 = {1.5, 1.5, 0.0, 0.0, 0};
    vec5d U2 = {0.5323, 0.3, 1.206, 0.0, 0};
    vec5d U3 = {0.138, 0.029, 1.206, 1.206, 0};
    vec5d U4 = {0.5323, 0.3, 0.0, 1.206, 0};
    vector<vector<vec5d>> ret;
    for (auto i = 0; i < N / 5; i++)
    {
        vector<vec5d> x_oriented_i(N);
        for (auto j = 0; j < 4 * N / 5; j++)
        {
            x_oriented_i[j] = U2;
        }
        for (auto j = 4 * N / 5; j < N; j++)
        {
            x_oriented_i[j] = U1;
        }
        ret.push_back(x_oriented_i);
    }
    for (auto i = N / 5; i < N; i++)
    {
        vector<vec5d> x_oriented_i(N);
        for (auto j = 0; j < 4 * N / 5; j++)
        {
            x_oriented_i[j] = U3;
        }
        for (auto j = 4 * N / 5; j < N; j++)
        {
            x_oriented_i[j] = U4;
        }
        ret.push_back(x_oriented_i);
    }
    return ret;
}

std::vector<std::vector<vec5d>> kernelLauncher::initialize_2Mach(int N)
{
    std::vector<std::vector<vec5d>> U(N, std::vector<vec5d>(4 * N, {1.4, 1, 0, 0, 0}));
    for (auto i = N - 1; i >= 0; i--)
    {
        double temp = double(N - i + 3) / sqrt(3);
        int l = N / 6 + int(temp);
        for (auto j = 0; j < l; j++)
        {
            U[i][j] = {8, 116.5, 4.125 * sqrt(3), -4.125, 0};
        }
    }
    return U;
}

std::vector<std::vector<vec5d>> kernelLauncher::initialize_jet800(int N)
{
    std::vector<std::vector<vec5d>> U(N, std::vector<vec5d>(3 * N, {0.14, 1, 0, 0, 0}));
    for (auto i = 0; i < N; ++i)
    {
        if (i > 9 * N / 10)
        {
            U[i][0] = {1.4, 1, 800, 0, 0};
        }
    }
    return U;
}

std::vector<std::vector<vec5d>> kernelLauncher::initialize_4cds(int N)
{
    vec5d U1 = {1.0, 1.0, 0.75, -0.5, 0};
    vec5d U2 = {2.0, 1.0, 0.75, 0.5, 0};
    vec5d U3 = {1.0, 1.0, -0.75, 0.5, 0};
    vec5d U4 = {3.0, 1.0, -0.75, -0.5, 0};
    vector<vector<vec5d>> ret;
    for (auto i = 0; i < N / 2; i++)
    {
        vector<vec5d> x_oriented_i(N);
        for (auto j = 0; j < N / 2; j++)
        {
            x_oriented_i[j] = U2;
        }
        for (auto j = N / 2; j < N; j++)
        {
            x_oriented_i[j] = U1;
        }
        ret.push_back(x_oriented_i);
    }
    for (auto i = N / 2; i < N; i++)
    {
        vector<vec5d> x_oriented_i(N);
        for (auto j = 0; j < N / 2; j++)
        {
            x_oriented_i[j] = U3;
        }
        for (auto j = N / 2; j < N; j++)
        {
            x_oriented_i[j] = U4;
        }
        ret.push_back(x_oriented_i);
    }
    return ret;
}
