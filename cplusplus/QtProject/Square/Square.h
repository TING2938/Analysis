#pragma once

#include <QtWidgets/QWidget>
#include "Grid.h"
#include <QTimer>
#include <QPainter>


class Square : public QWidget
{
	Q_OBJECT

public:
	Square(QWidget *parent = Q_NULLPTR);

	void keyReleaseEvent(QKeyEvent* event);
	
	void keyPressEvent(QKeyEvent* event);

	void paintEvent(QPaintEvent* event);
private:
	Grid* grid;
	QTimer* time;
};
