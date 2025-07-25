#pragma once
#include <QObject>
#include <fstream>
#include <QStringList>
#include <QImage>
#include "solvers.cuh"

#include <thread>
#include <QDebug>
#include <chrono>

class kernelLauncher : public QObject
{
Q_OBJECT // Qt's macro
    public : explicit kernelLauncher(QObject *parent = nullptr);
    QImage convertToQImage(const vector<vector<vec5d>> &data);
    static void writeToTxt(const QString &dir, const vector<vector<vec5d>> &data);

    double T = 0.0;
    int incrementCount = 50;
    QImage plot;
    QString outputDir = "./data";

    Q_INVOKABLE void launchWith(double CFL, double endT, double gamma, double alpha, int Nx, int Ny, int testID);
    Q_INVOKABLE void launchTest(double CFL, double endT, double gamma, double alpha, int Nx, int Ny, int testID)
    {
        for (int i = 0; i < incrementCount; i++)
        {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            emit progressChanged();
        }
        emit simulationFinished();
        std::this_thread::sleep_for(std::chrono::seconds(1));
        emit plotFinished(outputDir + "/data.png");
    }
    Q_INVOKABLE void reset();

signals:
    void simulationFinished();
    void plotFinished(const QString &filePath);
    void progressChanged();

private:
    std::vector<std::vector<vec5d>> initialize_4shocks(int N);
    std::vector<std::vector<vec5d>> initialize_2Mach(int N);
    std::vector<std::vector<vec5d>> initialize_jet800(int N);
    std::vector<std::vector<vec5d>> initialize_4cds(int N);

    const QVector<QRgb> RdBuColormap = {
        qRgb(103, 0, 31), qRgb(105, 0, 31), qRgb(108, 1, 31), qRgb(111, 2, 32),
        qRgb(114, 3, 32), qRgb(117, 4, 33), qRgb(120, 5, 33), qRgb(123, 6, 34),
        qRgb(126, 7, 34), qRgb(129, 8, 35), qRgb(132, 9, 35), qRgb(135, 10, 36),
        qRgb(138, 11, 36), qRgb(141, 12, 37), qRgb(144, 13, 37), qRgb(147, 14, 38),
        qRgb(150, 15, 38), qRgb(153, 16, 39), qRgb(155, 16, 39), qRgb(158, 17, 39),
        qRgb(161, 18, 40), qRgb(164, 19, 40), qRgb(167, 20, 41), qRgb(170, 21, 41),
        qRgb(173, 22, 42), qRgb(176, 23, 42), qRgb(178, 25, 43), qRgb(180, 28, 45),
        qRgb(181, 31, 46), qRgb(182, 33, 47), qRgb(184, 36, 49), qRgb(185, 39, 50),
        qRgb(187, 42, 51), qRgb(188, 45, 52), qRgb(190, 48, 54), qRgb(191, 50, 55),
        qRgb(192, 53, 56), qRgb(194, 56, 58), qRgb(195, 59, 59), qRgb(197, 62, 60),
        qRgb(198, 64, 62), qRgb(199, 67, 63), qRgb(201, 70, 65), qRgb(202, 73, 66),
        qRgb(204, 76, 67), qRgb(205, 79, 68), qRgb(206, 81, 70), qRgb(208, 84, 71),
        qRgb(209, 87, 73), qRgb(211, 90, 74), qRgb(212, 93, 75), qRgb(214, 96, 77),
        qRgb(215, 98, 79), qRgb(216, 101, 81), qRgb(217, 104, 83), qRgb(218, 106, 85),
        qRgb(219, 109, 87), qRgb(221, 112, 89), qRgb(222, 114, 91), qRgb(223, 117, 93),
        qRgb(224, 120, 95), qRgb(225, 123, 97), qRgb(226, 125, 99), qRgb(228, 128, 101),
        qRgb(229, 131, 104), qRgb(230, 133, 106), qRgb(231, 136, 108), qRgb(232, 139, 110),
        qRgb(234, 141, 112), qRgb(235, 144, 114), qRgb(236, 147, 116), qRgb(237, 150, 118),
        qRgb(238, 152, 120), qRgb(239, 155, 122), qRgb(241, 158, 124), qRgb(242, 160, 126),
        qRgb(243, 163, 128), qRgb(244, 166, 131), qRgb(244, 168, 134), qRgb(244, 170, 136),
        qRgb(245, 172, 139), qRgb(245, 174, 142), qRgb(245, 176, 144), qRgb(246, 178, 147),
        qRgb(246, 180, 150), qRgb(247, 182, 152), qRgb(247, 185, 155), qRgb(247, 187, 158),
        qRgb(248, 189, 161), qRgb(248, 191, 163), qRgb(248, 193, 166), qRgb(249, 195, 169),
        qRgb(249, 197, 171), qRgb(249, 199, 174), qRgb(250, 202, 177), qRgb(250, 204, 180),
        qRgb(250, 206, 182), qRgb(251, 208, 185), qRgb(251, 210, 188), qRgb(251, 212, 190),
        qRgb(252, 214, 193), qRgb(252, 216, 196), qRgb(253, 219, 199), qRgb(252, 220, 200),
        qRgb(252, 221, 202), qRgb(252, 222, 204), qRgb(252, 223, 206), qRgb(251, 224, 208),
        qRgb(251, 225, 210), qRgb(251, 226, 212), qRgb(251, 227, 214), qRgb(250, 228, 215),
        qRgb(250, 229, 217), qRgb(250, 231, 219), qRgb(250, 232, 221), qRgb(249, 233, 223),
        qRgb(249, 234, 225), qRgb(249, 235, 227), qRgb(249, 236, 229), qRgb(249, 237, 231),
        qRgb(248, 238, 232), qRgb(248, 239, 234), qRgb(248, 240, 236), qRgb(248, 242, 238),
        qRgb(247, 243, 240), qRgb(247, 244, 242), qRgb(247, 245, 244), qRgb(247, 246, 246),
        qRgb(246, 246, 246), qRgb(244, 245, 246), qRgb(243, 245, 246), qRgb(241, 244, 246),
        qRgb(240, 243, 245), qRgb(238, 243, 245), qRgb(237, 242, 245), qRgb(235, 241, 244),
        qRgb(234, 241, 244), qRgb(232, 240, 244), qRgb(231, 239, 244), qRgb(229, 238, 243),
        qRgb(228, 238, 243), qRgb(226, 237, 243), qRgb(225, 236, 243), qRgb(223, 236, 242),
        qRgb(222, 235, 242), qRgb(220, 234, 242), qRgb(219, 233, 241), qRgb(217, 233, 241),
        qRgb(216, 232, 241), qRgb(214, 231, 241), qRgb(213, 231, 240), qRgb(211, 230, 240),
        qRgb(210, 229, 240), qRgb(209, 229, 240), qRgb(206, 227, 239), qRgb(204, 226, 238),
        qRgb(201, 225, 237), qRgb(199, 223, 237), qRgb(196, 222, 236), qRgb(194, 221, 235),
        qRgb(191, 220, 235), qRgb(189, 218, 234), qRgb(186, 217, 233), qRgb(184, 216, 232),
        qRgb(181, 215, 232), qRgb(179, 213, 231), qRgb(176, 212, 230), qRgb(174, 211, 230),
        qRgb(171, 210, 229), qRgb(169, 208, 228), qRgb(167, 207, 228), qRgb(164, 206, 227),
        qRgb(162, 205, 226), qRgb(159, 203, 225), qRgb(157, 202, 225), qRgb(154, 201, 224),
        qRgb(152, 200, 223), qRgb(149, 198, 223), qRgb(147, 197, 222), qRgb(144, 196, 221),
        qRgb(141, 194, 220), qRgb(138, 192, 219), qRgb(135, 190, 218), qRgb(132, 188, 217),
        qRgb(128, 186, 216), qRgb(125, 184, 215), qRgb(122, 182, 214), qRgb(119, 180, 213),
        qRgb(116, 178, 211), qRgb(113, 176, 210), qRgb(110, 174, 209), qRgb(107, 172, 208),
        qRgb(104, 170, 207), qRgb(101, 168, 206), qRgb(97, 166, 205), qRgb(94, 164, 204),
        qRgb(91, 162, 203), qRgb(88, 160, 202), qRgb(85, 158, 201), qRgb(82, 156, 200),
        qRgb(79, 154, 199), qRgb(76, 152, 198), qRgb(73, 150, 197), qRgb(70, 148, 196),
        qRgb(67, 147, 195), qRgb(65, 145, 194), qRgb(64, 143, 193), qRgb(63, 141, 192),
        qRgb(61, 139, 191), qRgb(60, 138, 190), qRgb(59, 136, 189), qRgb(57, 134, 188),
        qRgb(56, 132, 187), qRgb(55, 131, 186), qRgb(53, 129, 185), qRgb(52, 127, 185),
        qRgb(51, 125, 184), qRgb(49, 124, 183), qRgb(48, 122, 182), qRgb(47, 120, 181),
        qRgb(45, 118, 180), qRgb(44, 117, 179), qRgb(43, 115, 178), qRgb(41, 113, 177),
        qRgb(40, 111, 176), qRgb(39, 109, 176), qRgb(37, 108, 175), qRgb(36, 106, 174),
        qRgb(35, 104, 173), qRgb(33, 102, 172), qRgb(32, 100, 170), qRgb(31, 98, 167),
        qRgb(30, 96, 164), qRgb(29, 94, 161), qRgb(28, 92, 158), qRgb(26, 90, 155),
        qRgb(25, 88, 152), qRgb(24, 86, 149), qRgb(23, 84, 147), qRgb(22, 81, 144),
        qRgb(21, 79, 141), qRgb(20, 77, 138), qRgb(19, 75, 135), qRgb(18, 73, 132),
        qRgb(17, 71, 129), qRgb(15, 69, 126), qRgb(14, 67, 123), qRgb(13, 64, 120),
        qRgb(12, 62, 117), qRgb(11, 60, 114), qRgb(10, 58, 111), qRgb(9, 56, 108),
        qRgb(8, 54, 105), qRgb(7, 52, 102), qRgb(6, 50, 99), qRgb(5, 48, 97)};
};