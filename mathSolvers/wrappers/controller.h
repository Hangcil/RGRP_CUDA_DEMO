#pragma once

#include <QObject>
#include <QThread>
#include "kernelLauncher.h"

class controller : public QObject
{
    Q_OBJECT
    Q_PROPERTY(double progress READ progress NOTIFY progressChanged)

public:
    explicit controller(QObject *parent = nullptr);
    ~controller();

    double progress() const;
    double m_progress = 0.0;

public slots:
    void startSimulation(double CFL, double endT, double gamma, double alpha, int Nx, int Ny, int testID);

private slots:
    void onProgressUpdated();
    void onSimulationFinished();
    void onPlotFinished(const QString &filePath);

signals:
    void operate(double CFL, double endT, double gamma, double alpha, int Nx, int Ny, int testID);
    void progressChanged();
    void simulationFinished();
    void plotFinished(const QString &filePath);

private:
    QThread m_workerThread;
    kernelLauncher *m_worker;
};
