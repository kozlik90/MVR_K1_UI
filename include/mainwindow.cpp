#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "testtask.h"
#include "maintask.h"
#include <qthread.h>
#include <QMessageBox>
#include "drawing.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow), testTask(nullptr), mainTask(nullptr), thread(nullptr)
{
    ui->setupUi(this);
    on_radioTest_btn_clicked();
    drawingTestReal = new Drawing(ui->widget1);
    drawingTestNumerical = new Drawing(ui->widget2);
    drawingTestInitApprox = new Drawing(ui->widget3);
    drawingTestError = new Drawing(ui->widget4);



    widget5 = new QWidget(ui->tab);
    widget6 = new QWidget(ui->tab_2);
    widget7 = new QWidget(ui->tab_3);
    widget8 = new QWidget(ui->tab_4);

    drawingMainNumerical1 = new Drawing(widget5);
    drawingMainNumerical2 = new Drawing(widget6);
    drawingMainInitApprox1 = new Drawing(widget7);
    drawingMainInitApprox2 = new Drawing(widget8);
    drawingMainError = new Drawing(ui->widget9);

    ui->widget1->raise();
    ui->widget2->raise();
    ui->widget3->raise();
    ui->widget4->raise();

}

MainWindow::~MainWindow()
{
    delete ui;
    delete testTask;
    delete thread;

}

void MainWindow::on_calcButton_clicked()

{
    delete testTask;
    testTask = nullptr;
    delete mainTask;
    mainTask = nullptr;
    delete thread;
    thread = nullptr;
    try{
        drawingTestReal->draw(testTask, Drawing::TestImageType::Clear);
        drawingTestNumerical->draw(testTask, Drawing::TestImageType::Clear);
        drawingTestInitApprox->draw(testTask, Drawing::TestImageType::Clear);
        drawingTestError->draw(testTask, Drawing::TestImageType::Clear);

        drawingMainNumerical1->draw(mainTask, Drawing::MainImageType::Clr);
        drawingMainNumerical2->draw(mainTask, Drawing::MainImageType::Clr);
        drawingMainInitApprox1->draw(mainTask, Drawing::MainImageType::Clr);
        drawingMainInitApprox2->draw(mainTask, Drawing::MainImageType::Clr);
        drawingMainError->draw(mainTask, Drawing::MainImageType::Clr);

        ui->calcButton->setEnabled(false);



        ui->infoBrowser->clear();
        ui->progressBar->setValue(0);
        thread = new QThread();

        int n = ui->n_Line->text().toInt();
        int m = ui->m_Line->text().toInt();
        double eps = ui->eps_Line->text().toDouble();
        double omega = ui->omega_Line->text().toDouble();
        int Nmax = ui->nMax_Line->text().toInt();
        int max = 0;


        if(ui->radioTest_btn->isChecked()) {
            ui->widget1->raise();
            ui->widget2->raise();
            ui->widget3->raise();
            ui->widget4->raise();

            testTask = new TestTask(n, m, eps, omega, Nmax);
            connect(thread, SIGNAL(started()), testTask, SLOT(compute()));
            max = 100 / eps;
            ui->progressBar->setRange(0, max);
            connect(testTask, SIGNAL(progress(int)), ui->progressBar, SLOT(setValue(int)));
            connect(testTask, SIGNAL(fileInfo(QString)), this, SLOT(print(QString)));
            connect(testTask, SIGNAL(info(QString)), this, SLOT(print(QString)));
            connect(testTask, SIGNAL(finished(bool)), ui->calcButton, SLOT(setEnabled(bool)));
            connect(thread, SIGNAL(finished()), this, SLOT(slotDraw()));
            testTask->moveToThread(thread);
        }
        else {
            widget5->raise();
            widget6->raise();
            widget7->raise();
            widget8->raise();

            mainTask = new MainTask(n, m, eps, omega, Nmax);
            connect(thread, SIGNAL(started()), mainTask, SLOT(compute()));
            max = 100 / (eps / 10.0);
            ui->progressBar->setRange(0, max);
            connect(mainTask, SIGNAL(progress(int)), ui->progressBar, SLOT(setValue(int)));
            connect(mainTask, SIGNAL(fileInfo(QString)), this, SLOT(print(QString)));
            connect(mainTask, SIGNAL(info(QString)), this, SLOT(print(QString)));
            connect(mainTask, SIGNAL(finished(bool)), ui->calcButton, SLOT(setEnabled(bool)));
            connect(mainTask, SIGNAL(finished(bool)), this, SLOT(slotDraw()), Qt::UniqueConnection);
            mainTask->moveToThread(thread);
        }
        thread->start();
        thread->quit();
    }
    catch(std::exception& ex) {
        QMessageBox box;
        box.setText(ex.what());
        box.setIcon(QMessageBox::Warning);
        box.exec();
        ui->calcButton->setEnabled(true);
    }



}

void MainWindow::slotDraw() {
    int cur = ui->tabWidget->currentIndex();
    ui->tabWidget->setCurrentIndex(0);
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setCurrentIndex(2);
    ui->tabWidget->setCurrentIndex(3);
    ui->tabWidget->setCurrentIndex(4);
    ui->tabWidget->setCurrentIndex(cur);

    widget5->setGeometry(ui->widget1->geometry());
    widget6->setGeometry(ui->widget2->geometry());
    widget7->setGeometry(ui->widget3->geometry());
    widget8->setGeometry(ui->widget4->geometry());

    drawingTestReal->resize(ui->widget1->size());
    drawingTestNumerical->resize(ui->widget2->size());
    drawingTestInitApprox->resize(ui->widget3->size());
    drawingTestError->resize(ui->widget4->size());

    drawingMainNumerical1->resize(widget5->size());
    drawingMainNumerical2->resize(widget6->size());
    drawingMainInitApprox1->resize(widget7->size());
    drawingMainInitApprox2->resize(widget8->size());
    drawingMainError->resize(ui->widget9->size());

    if(ui->radioTest_btn->isChecked()){
        drawingTestReal->draw(testTask, Drawing::TestImageType::Real);

        drawingTestNumerical->draw(testTask, Drawing::TestImageType::Numerical);

        drawingTestInitApprox->draw(testTask, Drawing::TestImageType::InitApprox);

        drawingTestError->draw(testTask, Drawing::TestImageType::Error);
    }
    else if (ui->radioMain_btn->isChecked()) {
        drawingMainNumerical1->draw(mainTask, Drawing::MainImageType::Numerical1);

        drawingMainNumerical2->draw(mainTask, Drawing::MainImageType::Numerical2);

        drawingMainInitApprox1->draw(mainTask, Drawing::MainImageType::InitApprox1);

        drawingMainInitApprox2->draw(mainTask, Drawing::MainImageType::InitApprox2);

        drawingMainError->draw(mainTask, Drawing::MainImageType::Err);
    }


}

void MainWindow::print(QString str) {
    ui->infoBrowser->append(str);
}


void MainWindow::on_radioTest_btn_clicked()
{

    ui->tabWidget->setTabText(0, "Истинное решение");
    ui->tabWidget->setTabText(1, "Численное решение");
    ui->tabWidget->setTabText(2, "Начальное приближение");
    ui->tabWidget->setTabText(3, "Погрешность");
    ui->tabWidget->setTabVisible(4, false);
    ui->tabWidget->setTabEnabled(4, 0);
}


void MainWindow::on_radioMain_btn_clicked()
{

    ui->tabWidget->setTabVisible(4, true);
    ui->tabWidget->setTabText(0, "Численное решение 1");
    ui->tabWidget->setTabText(1, "Численное решение 2");
    ui->tabWidget->setTabText(2, "Начальное приближение 1");
    ui->tabWidget->setTabText(3, "Начальное приближение 2");
    ui->tabWidget->setTabText(4, "Погрешность");
    ui->tabWidget->setTabEnabled(4, 1);
}

void MainWindow::resizeEvent(QResizeEvent *event)
{
    int cur = ui->tabWidget->currentIndex();
    ui->tabWidget->setCurrentIndex(0);
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setCurrentIndex(2);
    ui->tabWidget->setCurrentIndex(3);
    ui->tabWidget->setCurrentIndex(4);
    ui->tabWidget->setCurrentIndex(cur);

    drawingTestReal->resize(ui->widget1->size());
    drawingTestNumerical->resize(ui->widget2->size());
    drawingTestInitApprox->resize(ui->widget3->size());
    drawingTestError->resize(ui->widget4->size());

    widget5->setGeometry(ui->widget1->geometry());
    drawingMainNumerical1->resize(widget5->size());
    widget6->setGeometry(ui->widget2->geometry());
    drawingMainNumerical2->resize(widget6->size());
    widget7->setGeometry(ui->widget3->geometry());
    drawingMainInitApprox1->resize(widget7->size());
    widget8->setGeometry(ui->widget4->geometry());
    drawingMainInitApprox1->resize(widget8->size());
    drawingMainError->resize(ui->widget9->size());

    QMainWindow::resizeEvent(event);
}



