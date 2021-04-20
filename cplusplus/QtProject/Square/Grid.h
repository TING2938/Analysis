#pragma once

#include <QWidget>
#include <QLabel>
#include <QPoint>
#include <QKeyEvent>
#include <QPalette>

class Grid : public QWidget
{
	Q_OBJECT

public:
	
	Grid(QWidget *parent);
	~Grid();

	void setType(int t);

	void rotate(int dir);

	void moveLeft();
	void moveRight();

	void drop();

	void reset();

	void drawMap();

	void checkMap();

private:
	void setPos();

	int getRand(int beg, int end);

public:
	int a;
	int nrow;
	int ncol;
	int ntype;
	QPoint margin;
	int type;

private:
	QLabel* grid[4];
	QPoint pos[4];

	QVector<QPalette> palette;
	QVector<QVector<int>> map;
	QVector<QVector<QLabel*>> labelMap;
	QPoint typeVec[4];
};
