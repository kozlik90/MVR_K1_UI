#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qopenglwidget.h>

class TestTask;
class MainTask;
class Drawing;
QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_calcButton_clicked();

    void on_radioTest_btn_clicked();

    void on_radioMain_btn_clicked();


private:
    Ui::MainWindow *ui;
    TestTask* testTask;
    MainTask* mainTask;
    QThread* thread;

    Drawing* drawingTestReal;
    Drawing* drawingTestNumerical;
    Drawing* drawingTestInitApprox;
    Drawing* drawingTestError;

    Drawing* drawingMainNumerical1;
    Drawing* drawingMainNumerical2;
    Drawing* drawingMainInitApprox1;
    Drawing* drawingMainInitApprox2;
    Drawing* drawingMainError;

    QWidget* widget5;
    QWidget* widget6;
    QWidget* widget7;
    QWidget* widget8;
    void resizeEvent(QResizeEvent *event) override;
public slots:
    void print(QString);
    void slotDraw();
};
#endif // MAINWINDOW_H
