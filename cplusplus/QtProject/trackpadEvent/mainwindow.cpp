#include "mainwindow.h"
#include <QDebug>
#include <QEvent>
#include <iostream>
#include <QMouseEvent>
#include <QString>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), label(new QLabel)
{
    this->resize(1000, 800);
    label->setParent(this);
    label->resize(300, 100);
    label->move(300, 400);
    
}

MainWindow::~MainWindow()
{
}

void MainWindow::mouseMoveEvent(QMouseEvent* ev)
{
    auto lpos = ev->localPos();
    auto pos = ev->pos();
    auto gpos = ev->globalPos();
    QString str;
    str = QString("localPos: %1, %2\n").arg(lpos.x()).arg(lpos.y())
    + QString("pos: %1, %2\n").arg(pos.x()).arg(pos.y())
    + QString("globalPos: %1, %2\n").arg(gpos.x()).arg(gpos.y());
    label->setText(str);
}

