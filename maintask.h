#ifndef MAINTASK_H
#define MAINTASK_H

#include <QObject>
#include <exception>

class Drawing;

class MainTaskEx : public std::exception {
private:
    std::string msg;
public:
    MainTaskEx(std::string str) : msg(str) {}
    const char* what() const noexcept {return msg.c_str();}
};

class MainTask : public QObject
{
    Q_OBJECT
private:
    int n, n2, m, m2, iter, iter2;
    std::vector<std::vector<double>> v1, v2, Error;
    double eps, omega, maxDiff, maxDiff2, nevyazkaMax, nevyazkaMax2, maxError;
    double eps2, omega2;
    int Nmax, Nmax2;
    const double PI = 3.14159265358979323846;
    void set_GU(std::vector<std::vector<double>>&, int, int);
    void set_inter(std::vector<std::vector<double>>&, int, int);
    void calculateChisl(std::vector<std::vector<double>>&, int, int, double, double, int, int&, double&);
    void calculateError();
    double f_func(double, double);
    void writeChisl1();
    void writeChisl2();
    double time;
    friend class Drawing;
public:
    MainTask(int, int, double, double, int);
    void printInfo();
public slots:
    void compute();

signals:
    void progress(int);
    void fileInfo(QString);
    void info(QString);
    void finished(bool);
};

#endif // MAINTASK_H
