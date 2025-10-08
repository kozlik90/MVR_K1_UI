#include "drawing.h"
#include "maintask.h"
#include "qevent.h"
#include "testtask.h"


Drawing::Drawing(QWidget* pwgt) : QOpenGLWidget(pwgt), m_xRotate(0), m_yRotate(0), task(nullptr), mainTask(nullptr) {

}

void Drawing::draw(TestTask* t, TestImageType type)
{
    task = t;
    funcType = type;
    update();
}

void Drawing::draw(MainTask* t, MainImageType type) {
    mainTask = t;
    funcType = type;
    update();
}


void Drawing::initializeGL() {
    QOpenGLFunctions* pFunc = QOpenGLContext::currentContext()->functions();
    pFunc->glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glMatrixMode(GL_PROJECTION); // выбираем тип матрицы
    glLoadIdentity(); // загружаемся с этими настройками
    glOrtho(-18, 18, -18, 18, -160.0, 160.0); // устанавливем диапазон изменения координат
    glMatrixMode(GL_MODELVIEW);
}

void Drawing::resizeGL(int nWidth, int nHeight) {
    glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-1.0, 1.0, -1.0, 1.0, 1.0, 10.0);
}

void Drawing::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity(); // Сброс матрицы модели
    glScalef(m_scale, m_scale, m_scale);

    glTranslatef(m_xTrans, m_yTrans, 0.0f);
    glRotatef(35.264, 1.0, 0.0, 0.0);
    glRotatef(45, 0.0, 1.0, 0.0);
    glRotatef(-90, 1.0, 0.0, 0.0);


    glTranslatef(0.0f, 0.0f, -1.5f);

    glRotatef(m_yRotate, 0.0, 0.0, 1.0);
    glTranslatef(-0.5f, -0.5f, 0.0f);


    if(task != nullptr){
        switch (funcType) {
        case TestImageType::Real:
            drawAxes();
            drawTestReal();
            break;
        case TestImageType::Numerical:
            drawAxes();
            drawTestNumerical();
            break;
        case TestImageType::InitApprox:
            drawAxes();
            drawTestInitApprox();
            break;
        case TestImageType::Error:
            drawAxes();
            glTranslatef(0.0f, 0.0f, 1.0f);
            drawTestError();
            break;
        case TestImageType::Clear:
            clear();
            break;
        }


    }
    else if(mainTask != nullptr) {
        switch (funcType) {
        case MainImageType::Numerical1:
            drawAxes();
            glTranslatef(0.0f, 0.0f, 1.0f);
            drawMainNumerical1();
            break;
        case MainImageType::Numerical2:
            drawAxes();
            glTranslatef(0.0f, 0.0f, 1.0f);
            drawMainNumerical2();
            break;
        case MainImageType::InitApprox1:
            drawAxes();
            glTranslatef(0.0f, 0.0f, 1.0f);
            drawMainInitApprox1();
            break;
        case MainImageType::InitApprox2:
            drawAxes();
            glTranslatef(0.0f, 0.0f, 1.0f);
            drawMainInitApprox2();
            break;
        case MainImageType::Err:
            drawAxes();
            glTranslatef(0.0f, 0.0f, 1.0f);
            drawMainError();
            break;
        case MainImageType::Clr:
            clear();
            break;
        }
    }

    glFlush();
}

void Drawing::mousePressEvent(QMouseEvent* pe) {
    m_ptPosition = pe->pos();
}

void Drawing::mouseMoveEvent(QMouseEvent* pe) {
    if(pe->buttons() == Qt::RightButton) {

        m_xRotate += 180 * (GLfloat)(pe->y() - m_ptPosition.y()) / height();
        m_yRotate += 180 * (GLfloat)(pe->x() - m_ptPosition.x()) / width();
    }
    else {
        m_xTrans += (GLfloat)(pe->x() - m_ptPosition.x()) / width();
        m_yTrans += -(GLfloat)(pe->y() - m_ptPosition.y()) / height();
    }
    update();

    m_ptPosition = pe->pos();
}
void Drawing::wheelEvent(QWheelEvent* pe) {
    const qreal scaleCoef = 1.1;
    qreal newScale = pe->angleDelta().y() > 0 ? m_scale * scaleCoef : m_scale / scaleCoef;
    m_scale = newScale;
    update();
}


void Drawing::drawTestReal()
{
    glLineWidth(2.0f);
    int n = task->n;
    int m = task->m;
    std::vector<std::vector<double>> vec(n+1, std::vector<double>(m+1, 0));
    double hx = 1.0 / n;
    double ky = 1.0 / m;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;
            double z = vec[i][j] = task->u_func(x, y);
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (vec[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, vec[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (vec[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, vec[i][j]);
        }
        glEnd();
    }
    glFlush();
}

void Drawing::drawTestNumerical() {
    int n = task->n;
    int m = task->m;
    double hx = 1.0 / n;
    double ky = 1.0 / m;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            double z = task->u[i][j];
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (task->u[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, task->u[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (task->u[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, task->u[i][j]);
        }
        glEnd();
    }

    glFlush();
}

void Drawing::drawTestInitApprox() {
    int n = task->n;
    int m = task->m;
    TestTask test(n, m, task->eps, task->omega, task->Nmax);
    test.set_GU();
    test.set_inter();
    double hx = 1.0 / n;
    double ky = 1.0 / m;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            double z = test.u[i][j];
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (test.u[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, test.u[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (test.u[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, test.u[i][j]);
        }
        glEnd();
    }

    glFlush();
}

void Drawing::drawTestError() {
    int n = task->n;
    int m = task->m;
    double hx = 1.0 / n;
    double ky = 1.0 / m;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            double z = task->Error[i][j];
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (task->Error[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, task->Error[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (task->Error[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, task->Error[i][j]);
        }
        glEnd();
    }

    glFlush();
}

void Drawing::drawMainNumerical1() {
    int n = mainTask->n;
    int m = mainTask->m;
    double hx = 1.0 / n;
    double ky = 1.0 / m;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            double z = mainTask->v1[i][j];
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mainTask->v1[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mainTask->v1[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mainTask->v1[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mainTask->v1[i][j]);
        }
        glEnd();
    }

    glFlush();
}

void Drawing::drawMainNumerical2() {
    int n = mainTask->n;
    int m = mainTask->m;
    double hx = 1.0 / (2*n);
    double ky = 1.0 / (2*m);
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= 2*n; i++) {
        for (int j = 0; j <= 2*m; j++) {
            double z = mainTask->v2[i][j];
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= 2*n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= 2*m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mainTask->v2[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mainTask->v2[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= 2*m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= 2*n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mainTask->v2[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mainTask->v2[i][j]);
        }
        glEnd();
    }

    glFlush();
}

void Drawing::drawMainInitApprox1() {
    MainTask mt(mainTask->n, mainTask->m, mainTask->eps, mainTask->omega, mainTask->Nmax);
    mt.set_GU(mt.v1, mt.n, mt.m);
    mt.set_inter(mt.v1, mt.n, mt.m);
    int n = mt.n;
    int m = mt.m;
    double hx = 1.0 / n;
    double ky = 1.0 / m;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            double z = mt.v1[i][j];
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mt.v1[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mt.v1[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mt.v1[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mt.v1[i][j]);
        }
        glEnd();
    }

    glFlush();
}

void Drawing::drawMainInitApprox2() {
    MainTask mt(mainTask->n, mainTask->m, mainTask->eps, mainTask->omega, mainTask->Nmax);
    mt.set_GU(mt.v2, mt.n2, mt.m2);
    mt.set_inter(mt.v2, mt.n2, mt.m2);
    int n = mt.n2;
    int m = mt.m2;
    double hx = 1.0 / n;
    double ky = 1.0 / m;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            double z = mt.v2[i][j];
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mt.v2[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mt.v2[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mt.v2[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mt.v2[i][j]);
        }
        glEnd();
    }

    glFlush();
}

void Drawing::drawMainError() {
    int n = mainTask->n;
    int m = mainTask->m;
    double hx = 1.0 / n;
    double ky = 1.0 / m;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            double z = mainTask->Error[i][j];
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;

        }
    }

    for (int i = 0; i <= n; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j <= m; j++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mainTask->Error[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mainTask->Error[i][j]);

        }
        glEnd();
    }

    for (int j = 0; j <= m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= n; i++) {
            double x = i * hx;
            double y = j * ky;

            float t = (mainTask->Error[i][j] - minZ) / (maxZ - minZ);
            glColor3f(t, 0.0f, 1.0f - t);
            glVertex3d(x, y, mainTask->Error[i][j]);
        }
        glEnd();
    }

    glFlush();
}

void Drawing::clear()
{
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glFlush();
    update();
}


void Drawing::drawAxes(float length) {
    glLineWidth(2.0f); // можно увеличить толщину линий
    glBegin(GL_LINES);
    // Ось X (красный)
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(length, 0.0f, 1.0f);

    // Ось Y (зелёный)
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, length, 1.0f);

    // Ось Z (синий)
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, length + 1.0f);
    glEnd();
}
