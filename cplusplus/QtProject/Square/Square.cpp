#include "Square.h"
#include <QPushButton>
#include <QVector>

Square::Square(QWidget *parent)
	: QWidget(parent)
{
	grid = new Grid(this);
	
	this->resize(grid->a * (grid->ncol*2+1), grid->a * (grid->nrow + 10));
	time = new QTimer(this);
	time->setInterval(800);
	time->start();

	connect(time, &QTimer::timeout, grid, &Grid::drop);
}

void Square::keyReleaseEvent(QKeyEvent* event)
{
	switch (event->key())
	{
	case Qt::Key_Up:
		grid->rotate(1);
		break;

	case Qt::Key_Left:
		grid->moveLeft();
		break;

	case Qt::Key_Right:
		grid->moveRight();
		break;

	default:
		break;
	}
}

void Square::keyPressEvent(QKeyEvent* event)
{
	switch (event->key())
	{
	case Qt::Key_Down:
		grid->drop();
		break;

	default:
		break;
	}
}

void Square::paintEvent(QPaintEvent* event)
{
	QPainter painter(this);
	painter.setPen(QPen(QColor(0, 160, 230), 4));
	painter.setBrush(QColor(200, 200, 200));
	painter.drawRect(grid->margin.x(), grid->margin.y(), grid->ncol * grid->a, grid->nrow * grid->a);
}
