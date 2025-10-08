#ifndef DRAWING_H
#define DRAWING_H

#include <QOpenGLWidget>
#include "qopenglcontext.h"
#include "qopenglfunctions.h"


class TestTask;
class MainTask;
class Drawing : public QOpenGLWidget
{
private:
    GLuint m_nPyramid;
    GLfloat m_xRotate;
    GLfloat m_yRotate;
    QPoint m_ptPosition;
    GLfloat m_xTrans = 0;
    GLfloat m_yTrans = 0;
    GLfloat m_scale = 1;

    TestTask* task;
    MainTask* mainTask;

    const double PI = 3.14159265358979323846;
    int funcType;

protected:
    virtual void initializeGL();
    virtual void resizeGL(int nWidth, int nHeight);
    virtual void paintGL();
    virtual void mousePressEvent(QMouseEvent* pe);
    virtual void mouseMoveEvent(QMouseEvent* pe);
    virtual void wheelEvent(QWheelEvent* pe);
    void drawAxes(float length = 1.5f);

public:
    Drawing(QWidget* pwgt = 0);
    enum TestImageType {Real = 1, Numerical, InitApprox, Error, Clear};
    enum MainImageType {Numerical1 = 1, Numerical2, InitApprox1, InitApprox2, Err, Clr};
    void draw(TestTask*, TestImageType);
    void draw(MainTask*, MainImageType);
    void drawTestReal();
    void drawTestNumerical();
    void drawTestInitApprox();
    void drawTestError();

    void drawMainNumerical1();
    void drawMainNumerical2();
    void drawMainInitApprox1();
    void drawMainInitApprox2();
    void drawMainError();

    void clear();
};

#endif // DRAWING_H
