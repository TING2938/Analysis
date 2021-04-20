#pragma once

#include <QtWidgets/QWidget>
#include "Grid.h"

#define SIZE 5

class Game : public QWidget
{
	Q_OBJECT

public:
	Game(QWidget *parent = Q_NULLPTR);

	void addDefaultGrid();

	void merge(int direction);

private:
	int getRand(int beg, int end);

	void mergeFunc(int row, bool col);
	void mergeFunc(bool row, int col);

	void keyReleaseEvent(QKeyEvent* event);

private:
	Grid* grid[SIZE][SIZE];
};
