#include "testtask.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <iomanip>
#include <QtCore>
#include <ctime>

using namespace std;

TestTask::TestTask(int n_, int m_, double eps_, double omega_, int Nmax_) :  n(n_),
    m(m_), u(n+1, std::vector<double>(m+1, 0)), Error(n+1, std::vector<double>(m+1, 0)), eps(eps_), omega(omega_), Nmax(Nmax_), hx(1.0/n), ky(1.0/m),
    iter(0), maxDiff(0.0), nevyazkaMax(0.0), maxError(0.0), time(0.0)
{
    if(n <= 0)
        throw TestTaskEx("n must be >= 0");
    else if (m <= 0)
        TestTaskEx("m must be >= 0");
    else if(eps < 0)
        throw TestTaskEx("eps must be > 0");
    else if (omega < 0 || omega > 2)
        throw TestTaskEx("omega must be [0, 2), if you choose omega = 2, then omega will take the optimal value for the given n and m");
    else if(Nmax < 0)
        throw TestTaskEx("Nmax must be > 0");
    if (omega >= 2) {
        omega = 2.0 / (1 + 2 * sin(M_PI * hx / 2));
    }
}



double TestTask::u_func(const double x, const double y) {
    return exp(pow(sin(M_PI * x * y), 2));
}

double TestTask::f_func(const double x, const double y) {
    return -2 * M_PI * M_PI * exp(pow(sin(M_PI * x * y), 2)) * (y * y * pow(cos(M_PI * x * y), 2) + x * x * pow(cos(M_PI * x * y), 2) + pow(sin(M_PI * x * y), 2) * cos(2 * M_PI * x * y) * (y * y + x * x));
}

void TestTask::set_GU() {

    for (int i = 0; i <= n; i++) {
        double x = i * hx;
        u[i][0] = u_func(x, 0.0);
        u[i][m] = u_func(x, 1.0);
    }

    for (int j = 0; j <= m; j++) {
        double y = j * ky;
        u[0][j] = u_func(0.0, y);
        u[n][j] = u_func(1.0, y);
    }
}

void TestTask::set_inter() {
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            double x = i * hx;
            // linear interpolation x:
            u[i][j] = u[0][j] * (1 - x) + u[n][j] * x; //
        }
    }
}

void TestTask::calculateIst() {
    QString str = "";
    std::ofstream outFile1("1istinnoe.csv");
    if (!outFile1.is_open()) {
        str = "File 1 is not open.\n";
    }
    else {
    outFile1 << ";";
    for (int i = 0; i <= n; i++) {
        outFile1 << "x" << i;
        if (i < n)
            outFile1 << ";";
    }
    outFile1 << "\n";

    for (int j = 0; j < m + 1; j++) {
        outFile1 << "y" << j << ";";
        for (int i = 0; i < n + 1; i++) {
            double x = j * hx;
            double y = i * ky;
            outFile1 << std::setprecision(30) << u_func(x, y);
            if (i < n)
                outFile1 << ";";

        }
        outFile1 << "\n";
    }
    outFile1.close();
    str = "File 1istinnoe.csv is ready\n";
    }
    emit fileInfo(str);
}

void TestTask::calculateChisl() {
    double pow_h2 = 1.0 / (hx * hx);
    double pow_k2 = 1.0 / (ky * ky);
    double denom = 2 * pow_h2 + 2 * pow_k2;
    do {
        maxDiff = 0.0;

        for (int i = 1; i < n; i++) {
            for (int j = 1; j < m; j++) {

                double x = i * hx;
                double y = j * ky;

                double rhs = f_func(x, y) + pow_h2 * (u[i - 1][j] + u[i + 1][j]) + pow_k2 * (u[i][j - 1] + u[i][j + 1]);
                double u_old = u[i][j];
                double u_new = (1 - omega) * u_old + (omega / denom) * rhs;
                u[i][j] = u_new;
                double diff = fabs(u_new - u_old);
                if (diff > maxDiff) maxDiff = diff;
            }
        }

        iter++;
        emit progress(100 / maxDiff);
    } while (maxDiff > eps && iter < Nmax);

    std::ofstream outFile2("2chislennoe.csv");
    QString str = "";
    if (!outFile2.is_open()) {
        str = "File 2 is not open.\n";
    }
    else {
        outFile2 << ";";
        for (int i = 0; i <= n; i++) {
            outFile2 << "x" << i;
            if (i < n)
                outFile2 << ";";
        }
        outFile2 << "\n";

        for (int j = 0; j < m + 1; j++) {
            outFile2 << "y" << j << ";";
            for (int i = 0; i < n + 1; i++) {
                outFile2 << std::setprecision(30) << u[i][j];
                if (i < n)
                    outFile2 << ";";

            }
            outFile2 << "\n";
        }
        outFile2.close();
        str = "File 2chislennoe.csv is ready\n";
    }
    emit fileInfo(str);

}

void TestTask::calculateError() {
    double pow_h2 = 1.0 / (hx * hx);
    double pow_k2 = 1.0 / (ky * ky);

    for (int j = 0; j < m + 1; j++) {
        double y = j * ky;
        for (int i = 0; i < n + 1; i++) {
            double x = i * hx;
            Error[i][j] = fabs(u[i][j] - u_func(x, y));
            if (Error[i][j] > maxError) {
                maxError = Error[i][j];
            }
        }
    }
    QString str = "";
    std::ofstream outFile3("3pogreshnost.csv");
    if (!outFile3.is_open()) {
        str = "File 3 is not open.\n";

    }
    else {
    outFile3 << ";";
    for (int i = 0; i <= n; i++) {
        outFile3 << "x" << i;
        if (i < n)
            outFile3 << ";";
    }
    outFile3 << "\n";

    for (int j = 0; j < m + 1; j++) {
        outFile3 << "y" << j << ";";
        for (int i = 0; i < n + 1; i++) {
            outFile3 << std::setprecision(30) << Error[i][j];
            if (i < n)
                outFile3 << ";";
        }
        outFile3 << "\n";
    }

    outFile3.close();
    str = "File 3pogreshnost.csv is ready\n";
    }
    emit fileInfo(str);

    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            double x = i * hx;
            double y = j * ky;
            double r = -pow_h2 * u[i - 1][j] - pow_k2 * u[i][j - 1]
                       + (2 * pow_h2 + 2 * pow_k2) * u[i][j]
                       - pow_h2 * u[i + 1][j] - pow_k2 * u[i][j + 1]
                       - f_func(x, y);
            if (fabs(r) > nevyazkaMax) {
                nevyazkaMax = fabs(r);
            }
        }
    }
    emit progress(100 / eps);
}
void TestTask::compute() {
    unsigned int start = clock();
    set_GU();
    set_inter();
    calculateIst();
    calculateChisl();
    calculateError();
    unsigned int end = clock();
    time = (end - start) / 1000.;
    printInfo();
    emit finished(true);

}

void TestTask::printInfo() {
    QString str = "";
    str += "The number of partitions by x: n = " + QString::number(n) + " and number of partitions by y: m = " + QString::number(m) + "\n";
    str += "The number of partitions by x: n = " + QString::number(n) + " and number of partitions by y: m = " + QString::number(m) + "\n";
    str += "SOR with the parameter omega = " + QString::number(omega) + ", criteria for stopping by accuracy eps = " + QString::number(eps) + " and by the number of iterations Nmax = " + QString::number(Nmax) + "\n";
    str += "iterations spent N = " + QString::number(iter) + "; the achieved accuracy of the iterative method eps^(N) = " + QString::number(maxDiff) + "\n";
    str += "The discrepancy ||R^(N)|| = " + QString::number(nevyazkaMax) + "\n";
    str += "The test task should be solved with an error of no more than eps = 0.5*10^(-6); the task was solved with an error of eps1 = " + QString::number(maxError) + "\n";
    str += "x interpolation is used as an initial approximation.\n";
    str += QString::number(time) + "\n";
    emit info(str);
}

