#include "mainwindow.h"
#include <QChartView>
#include <QChart>
#include <QLineSeries>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFileDialog>
#include <cmath>
#include <itp/fileio>

using namespace QtCharts;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    this->resize(1000, 800);
    line = new QLineSeries();
    chart = new QChart();
    cv = new QChartView();
    this->setUI();
    this->initBtn();
}

MainWindow::~MainWindow()
{
}

void MainWindow::setUI()
{
    QWidget* mainWidget = new QWidget;
    this->setCentralWidget(mainWidget);

    QHBoxLayout* Frame = new QHBoxLayout();
    mainWidget->setLayout(Frame);

    cv = new QChartView();
    cv->resize(600, 600);
    Frame->addWidget(cv);
    
    QWidget* right = new QWidget();
    Frame->addWidget(right);
    QVBoxLayout* rightLayout = new QVBoxLayout();
    right->setLayout(rightLayout);
    btn_open = new QPushButton(tr("Open"));
    btn_plot = new QPushButton(tr("Plot"));
    btn_save = new QPushButton(tr("Save"));
    btn_clear = new QPushButton(tr("Clear"));
    rightLayout->addWidget(btn_open);
    rightLayout->addWidget(btn_plot);
    rightLayout->addWidget(btn_save);
    rightLayout->addWidget(btn_clear);
}

void MainWindow::initBtn()
{
    connect(btn_open, &QPushButton::clicked, [=]{
        fnm = QFileDialog::getOpenFileName(this, "Open");
        if (!fnm.isEmpty())
        {
            this->loadData();
            this->plot();
        }
    });

    connect(btn_save, &QPushButton::clicked, [=]{
        // save
        auto str = QFileDialog::getSaveFileName(this, tr("Save picture"), "",
                                                tr("PNG file (*.png) ;; JPG file (*.jpg)"));
        cv->grab().save(str, "PNG");
    });
}

void MainWindow::loadData()
{
    auto tmp = itp::loadtxt(fnm.toStdString(), Eigen::Dynamic, 2, "#@");
    for (size_t i = 0; i != tmp.rows(); ++i)
    {
        line->append(tmp(i, 0), tmp(i, 1));
    }
}

void MainWindow::plot()
{
    chart->addSeries(line);
    chart->createDefaultAxes();
    chart->legend()->hide();

    cv->setChart(chart);

}
