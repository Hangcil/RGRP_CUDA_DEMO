#include "controller.h"

controller::controller(QObject *parent) : QObject(parent)
{
    m_worker = new kernelLauncher;
    m_worker->moveToThread(&m_workerThread);

    // Connect signals for thread management and communication
    connect(&m_workerThread, &QThread::finished, m_worker, &QObject::deleteLater);
    connect(this, &controller::operate, m_worker, &kernelLauncher::launchWith);
    connect(m_worker, &kernelLauncher::progressChanged, this, &controller::onProgressUpdated);
    connect(m_worker, &kernelLauncher::simulationFinished, this, &controller::onSimulationFinished);
    connect(m_worker, &kernelLauncher::plotFinished, this, &controller::onPlotFinished);

    m_workerThread.start();
}

controller::~controller()
{
    m_workerThread.quit();
    m_workerThread.wait();
}

double controller::progress() const
{
    return m_progress;
}

void controller::startSimulation(double CFL, double endT, double gamma, double alpha, int Nx, int Ny, int testID)
{
    m_progress = 0.0;
    emit progressChanged();
    emit operate(CFL, endT, gamma, alpha, Nx, Ny, testID);
}

void controller::onProgressUpdated()
{
    m_progress += 1.0 / double(m_worker->incrementCount);
    emit progressChanged();
}

void controller::onSimulationFinished()
{
    emit simulationFinished();
}

void controller::onPlotFinished(const QString &filePath)
{
    emit plotFinished(filePath);
}
