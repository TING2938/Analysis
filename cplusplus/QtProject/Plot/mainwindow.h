#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLineSeries>
#include <QPushButton>
#include <QChart>
#include <QChartView>
#include <QLabel>
#include <QLineEdit>
#include <QString>

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void setUI();
    void initBtn();
    void loadData();
    void plot();

private:
    QString fnm;
    QtCharts::QLineSeries* line;
    QtCharts::QChart* chart;
    QtCharts::QChartView* cv;
    QPushButton* btn_plot;
    QPushButton* btn_clear;
    QPushButton* btn_open;
    QPushButton* btn_save;

};
#endif // MAINWINDOW_H
