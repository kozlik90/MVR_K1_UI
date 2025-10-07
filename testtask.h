#ifndef TESTTASK_H
#define TESTTASK_H

#include <vector>
#include <QObject>
#include <exception>

class Drawing;

class TestTaskEx : public std::exception {
private:
    std::string msg;
public:
    TestTaskEx(std::string str) : msg(str) {}
    const char* what() const noexcept override {return msg.c_str();}
};

class TestTask : public QObject
{
    Q_OBJECT
private:
    int n, m, iter;
    double hx, ky, eps, omega, maxDiff, nevyazkaMax, maxError;
    int Nmax;
    double time;
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> Error;
    friend class Drawing;

public:
    TestTask(int, int, double, double, int);
    void printInfo();

private:
    static double u_func(double x, double y);
    double f_func(double x, double y);
    void set_GU();
    void set_inter();
    void calculateIst();
    void calculateChisl();
    void calculateError();

public slots:
    void compute();

signals:
    void progress(int);
    void fileInfo(QString);
    void info(QString);
    void finished(bool);

};

#endif // TESTTASK_H
