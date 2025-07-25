#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QLibraryInfo>
#include <QQmlContext>
#include <QThread>
#include "../mathSolvers/wrappers/controller.h"

int main(int argc, char *argv[])
{
    QGuiApplication app(argc, argv);

    QQmlApplicationEngine engine;

    auto launcher = new controller();
    engine.rootContext()->setContextProperty("launcher", launcher);

    QObject::connect(
        &engine,
        &QQmlApplicationEngine::objectCreationFailed,
        &app,
        []()
        { QCoreApplication::exit(-1); },
        Qt::QueuedConnection);
    engine.loadFromModule("RGRP_QML", "Main");

    return app.exec();
}
