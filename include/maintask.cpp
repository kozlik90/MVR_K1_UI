#include "maintask.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <qdebug.h>

MainTask::MainTask(int n_, int m_, double eps_, double omega_, int Nmax_) : n(n_), m(m_), n2(2*n_), m2(2*m_), v1(n+1, std::vector<double>(m+1, 0)),
    v2(n2+1, std::vector<double>(m2 + 1, 0)), Error(n+1, std::vector<double>(m+1, 0)), eps(eps_), eps2(eps_/10.0), omega(omega_),
    omega2(omega_), Nmax(Nmax_), Nmax2(Nmax_ * 1.5), iter(0), iter2(0), maxDiff(0.0), maxDiff2(0.0), nevyazkaMax(0.0), nevyazkaMax2(0.0), time(0.0)
{
    if(n <= 0)
        throw MainTaskEx("n must be >= 0");
    else if (m <= 0)
        MainTaskEx("m must be >= 0");
    else if(eps < 0)
        throw MainTaskEx("eps must be > 0");
    else if(omega < 0 || omega > 2)
        throw MainTaskEx("omega must be [0, 2), if you choose omega = 2, then omega will take the optimal value for the given n and m");
    else if(Nmax < 0)
        throw MainTaskEx("Nmax must be > 0");
    if(omega >= 2) {
        double hx = 1.0 / n;
        double hx2 = 1.0 / n2;
        omega = 2.0 / (1 + 2 * sin(PI * hx / 2));
        omega2 = 2.0 / (1 + 2 * sin(PI * hx2 / 2));
    }
}

double MainTask::f_func(const double x, const double y) {
    return sin(PI * x * y) * sin(PI * x * y);
}

void MainTask::set_GU(std::vector<std::vector<double>>& u,  const int n_, const int m_) {
    double hx = 1.0 / n_;
    double ky = 1.0 / m_;
    for (int j = 0; j <= m_; j++) {
        double y = j * ky;
        u[0][j] = sin(PI * y);  // x = 0
        u[n_][j] = sin(PI * y);  // x = 1
    }

    for (int i = 0; i <= n_; i++) {
        double x = i * hx;
        u[i][0] = x - x * x;      // y = 0
        u[i][m_] = x - x * x;      // y = 1
    }
}

void MainTask::set_inter(std::vector<std::vector<double>>& u, const int n_, const int m_) {
    double hx = 1.0 / n_;
    for (int i = 1; i < n_; i++) {
        double x = i * hx;
        for (int j = 1; j < m_; j++) {
            u[i][j] = u[0][j] * (1 - x) + u[n_][j] * x;
        }
    }
}

void MainTask::calculateChisl(std::vector<std::vector<double>>& v, const int n_, const int m_, const double eps_, const double omega_, const int Nmax_, int& iter_, double& maxDiff_) {
    double hx = 1.0 / n_;
    double ky = 1.0 / m_;

    double pow_h2 = 1.0 / (hx * hx);
    double pow_k2 = 1.0 / (ky * ky);
    double denom = 2 * pow_h2 + 2 * pow_k2;

    iter_ = 0;
    maxDiff_ = 0.0;

    do {
        maxDiff_ = 0.0;

        for (int i = 1; i < n_; i++) {
            for (int j = 1; j < m_; j++) {
                double x = i * hx;
                double y = j * ky;

                double rhs = f_func(x, y) + pow_h2 * (v[i - 1][j] + v[i + 1][j]) + pow_k2 * (v[i][j - 1] + v[i][j + 1]);
                double v_old = v[i][j];
                double v_new = (1 - omega_) * v_old + (omega_ / denom) * rhs;
                v[i][j] = v_new;
                double diff = fabs(v_new - v_old);
                if (diff > maxDiff_) {
                    maxDiff_ = diff;
                }
            }
        }
        iter_++;
        if(n_ == n2) {
            emit progress(100 / maxDiff_);
        }

    } while (maxDiff_ > eps_ && iter_ < Nmax_);

    if(n_ == n) {
        writeChisl1();
    }
    else if (n_ == n2) writeChisl2();
}

void MainTask::calculateError() {
    double hx2 = 1.0 / n;
    double ky2 = 1.0 / m;
    double pow2_h2 = 1.0 / (hx2 * hx2);
    double pow2_k2 = 1.0 / (ky2 * ky2);
    double h2x = 1.0 / n2;
    double k2y = 1.0 / m2;

    for (int j = 0; j < m + 1; j++) {
        for (int i = 0; i < n + 1; i++) {
            Error[i][j] = fabs(v1[i][j] - v2[2 * i][2 * j]);
            if (Error[i][j] > maxError) {
                maxError = Error[i][j];
            }
        }
    }

    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            double x = i * hx2;
            double y = j * ky2;
            double r = -pow2_h2 * v1[i - 1][j] - pow2_k2 * v1[i][j - 1]
                       + (2 * pow2_h2 + 2 * pow2_k2) * v1[i][j]
                       - pow2_h2 * v1[i + 1][j] - pow2_k2 * v1[i][j + 1]
                       - f_func(x, y);
            if (fabs(r) > nevyazkaMax)
                nevyazkaMax = fabs(r);
        }
    }
    double powh2 = 1.0 / (h2x * h2x);
    double powk2 = 1.0 / (k2y * k2y);
    nevyazkaMax2 = 0.0;
    for (int i = 1; i < n2; i++) {
        for (int j = 1; j < m2; j++) {
            double x = i * h2x;
            double y = j * k2y;
            double r = -powh2 * v2[i - 1][j] - powk2 * v2[i][j - 1]
                       + (2 * powh2 + 2 * powk2) * v2[i][j]
                       - powh2 * v2[i + 1][j] - powk2 * v2[i][j + 1]
                       - f_func(x, y);
            if (fabs(r) > nevyazkaMax2)
                nevyazkaMax2 = fabs(r);
        }
    }
    writeError();
    emit progress(100 / eps2);
}

void MainTask::writeChisl1() {
    QString str = "";
    std::ofstream outFile11("1_chislennoe.csv");
    if (!outFile11.is_open()) {
        str = "File 1_chislennoe is not open.\n";
    }
    else {
        outFile11 << ";";
        for (int i = 0; i <= n; i++) {
            outFile11 << "x" << i;
            if (i < n)
                outFile11 << ";";
        }
        outFile11 << "\n";
        for (int j = 0; j < m + 1; j++) {
            outFile11 << "y" << j << ";";
            for (int i = 0; i < n + 1; i++) {
                outFile11 << std::setprecision(30) << v1[i][j];
                if (i < n)
                    outFile11 << ";";
            }
            outFile11 << "\n";
        }
        outFile11.close();
        str = "File 1_chislennoe.csv is ready\n";
    }
    emit fileInfo(str);
}

void MainTask::writeChisl2() {
    QString str = "";
    std::ofstream outFile22("2_chislennoe.csv");
    if (!outFile22.is_open()) {
        str = "File 2_chislennoe is not open.\n";
    }
    else {
        outFile22 << ";";
        for (int i = 0; i <= n2; i++) {
            outFile22 << "x" << i*0.5;
            if (i < n2)
                outFile22 << ";";
        }
        outFile22 << "\n";
        for (int j = 0; j < m2 + 1; j++) {
            outFile22 << "y" << j*0.5 << ";";
            for (int i = 0; i < n2 + 1; i++) {
                outFile22 << std::setprecision(30) << v2[i][j];
                if (i < n2)
                    outFile22 << ";";

            }
            outFile22 << "\n";
        }
        outFile22.close();
        str = "File 2_chislennoe2.csv is ready\n";
    }
    emit fileInfo(str);
}

void MainTask::writeError() {
    QString str = "";
    std::ofstream outFile33("3_pogreshnost.csv");
    if (!outFile33.is_open()) {
        str = "File 3_pogreshnost is not open.\n";
    }
    else {
        outFile33 << ";";
        for (int i = 0; i <= n; i++) {
            outFile33 << "x" << i;
            if (i < n)
                outFile33 << ";";
        }
        outFile33 << "\n";

        for (int j = 0; j < m + 1; j++) {
            outFile33 << "y" << j << ";";
            for (int i = 0; i < n + 1; i++) {

                outFile33 << std::setprecision(30) << Error[i][j];
                if (i < n)
                    outFile33 << ";";
            }
            outFile33 << "\n";
        }
        outFile33.close();
        str = "File 3_pogreshnost.csv is ready\n";
    }
    emit fileInfo(str);
}





void MainTask::compute() {
    unsigned int start = clock();
    set_GU(v1, n, m);
    set_inter(v1, n, m);
    set_GU(v2, n2, m2);
    set_inter(v2, n2, m2);
    calculateChisl(v1, n, m, eps, omega, Nmax, iter, maxDiff);
    calculateChisl(v2, n2, m2, eps2, omega2, Nmax2, iter2, maxDiff2);
    calculateError();
    unsigned int end = clock();
    time = (end - start) / 1000.;
    printInfo();
    emit finished(true);
}

void MainTask::printInfo() {
    QString str = "";
    str += "The number of partitions by x: n = " + QString::number(n) + " and number of partitions by y: m = " + QString::number(m) + "\n";
    str += "SOR with the parameter omega = " + QString::number(omega) + ", criteria for stopping by accuracy eps = " + QString::number(eps) + " and by the number of iterations Nmax = " + QString::number(Nmax) + "\n";
    str += "iterations spent N = " + QString::number(iter) + " the achieved accuracy of the iterative method eps^(N) = " + QString::number(maxDiff) + "\n";
    str += "The discrepancy ||R^(N)|| =" + QString::number(nevyazkaMax) + "\n";
    str += "A half-step grid is used to control the accuracy, SOR with the parameter omega2 = " + QString::number(omega2) + ", criteria for stopping by accuracy eps2 = " + QString::number(eps2) + " and by the number of iterations Nmax2 " + QString::number(Nmax2) + "\n";
    str += "iterations spent  N2 = " + QString::number(iter2) + " the achieved accuracy of the iterative method eps^(N2) = " + QString::number(maxDiff2) + "\n";
    str += "The discrepancy ||R^(N2)|| = " + QString::number(nevyazkaMax2) + "\n";
    str += "The main task should be solved with an error of no more than eps = 0.5*10^(-6); the task was solved with an error of eps2 = " + QString::number(maxError) + "\n";
    str += "x interpolation is used as an initial approximation.\n";
    str += QString::number(time);
    emit info(str);
}
