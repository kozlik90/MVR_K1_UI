QT       += core gui

QT += openglwidgets
QT += opengl
QT += widgets
QT += concurrent

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

windows {
LIBS += -lopengl32
}

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    include/main.cpp \
    include/testtask.cpp \
    include/maintask.cpp \
    include/mainwindow.cpp \
    include/drawing.cpp


HEADERS += \
    include/testtask.h \
    include/maintask.h \
    include/mainwindow.h \
    include/drawing.h \

FORMS += \
    include/mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
