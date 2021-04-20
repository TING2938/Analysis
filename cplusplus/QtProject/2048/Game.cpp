#include "Game.h"
#include <QVector>
#include <QList>
#include <QString>
#include <QPushButton>
#include <QKeyEvent>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QSpacerItem>
#include <random>
#include <algorithm>

Game::Game(QWidget *parent)
	: QWidget(parent)
{
	this->setFixedSize(SIZE*120+420, SIZE*120+120);
	this->setWindowTitle("2048");
	
	auto leftLayout = new QVBoxLayout();
	QGridLayout* gridLayout = new QGridLayout();
	leftLayout->addStretch(1);
	leftLayout->addLayout(gridLayout);
	leftLayout->addStretch(1);

	gridLayout->setSpacing(6);
	for (int i = 0; i != SIZE; ++i)
	{
		for (int j = 0; j != SIZE; ++j)
		{
			grid[i][j] = new Grid();
			gridLayout->addWidget(grid[i][j], j, i);
			grid[i][j]->setValue(0);
			grid[i][j]->draw();
		}
	}

	auto rightLayout = new QVBoxLayout();
	auto btn_up = new QPushButton("up");
	auto btn_down = new QPushButton("down");
	auto btn_left = new QPushButton("left");
	auto btn_right = new QPushButton("right");
	auto btn_reset = new QPushButton("reset");

	rightLayout->addStretch(2);
	rightLayout->addWidget(btn_up);
	auto hbox = new QHBoxLayout();
	hbox->addWidget(btn_left);
	hbox->addWidget(btn_right);
	rightLayout->addLayout(hbox);
	rightLayout->addWidget(btn_down);
	rightLayout->addStretch(3);
	rightLayout->addWidget(btn_reset);
	rightLayout->addStretch(2);
	
	auto layout = new QHBoxLayout();
	layout->addStretch(2);
	layout->addLayout(leftLayout);
	layout->addStretch(3);
	layout->addLayout(rightLayout);
	layout->addStretch(3);

	this->setLayout(layout);

	connect(btn_reset, &QPushButton::clicked, [=]() {
		for (int i = 0; i != SIZE; ++i)
		{
			for (int j = 0; j != SIZE; ++j)
			{
				grid[i][j]->setValue(0);
				grid[i][j]->draw();
			}
		}
	});

	connect(btn_up, &QPushButton::clicked, [=]() {
		this->merge(0);
	});
	connect(btn_down, &QPushButton::clicked, [=]() {
		this->merge(1);
	});
	connect(btn_left, &QPushButton::clicked, [=]() {
		this->merge(2);
	});
	connect(btn_right, &QPushButton::clicked, [=]() {
		this->merge(3);
	});
	
}

void Game::addDefaultGrid()
{
	QList<QPoint> zeroValue;
	for (int i = 0; i != SIZE; ++i)
	{
		for (int j = 0; j != SIZE; ++j)
		{
			if (grid[i][j]->getValue() == 0)
			{
				zeroValue.append(QPoint(i, j));
			}
		}
	}
	
	if (!zeroValue.empty())
	{
		auto point = zeroValue[getRand(0, zeroValue.size()-1)];
		grid[point.x()][point.y()]->setValue(getRand(1, 2)*2);
		grid[point.x()][point.y()]->draw();
	}
	
}

void Game::merge(int direction)
{
	bool up_down = (direction < 2);
	bool reverse = (direction % 2);

	for (int i = 0; i < SIZE; ++i)
	{
		if (up_down)
		{
			mergeFunc(i, reverse);
		}
		else
		{
			mergeFunc(reverse, i);
		}
	}

	this->addDefaultGrid();
}

int Game::getRand(int beg, int end)
{
	static std::random_device rd;
	static std::default_random_engine e(rd());
	std::uniform_int_distribution<int> dist(beg, end);
	return dist(e);
}

void Game::mergeFunc(int row, bool col)
{
	QList<int> temp;
	for (int i = 0; i != SIZE; ++i)
	{
		if (grid[row][i]->getValue() != 0)
		{
			temp.append(grid[row][i]->getValue());
		}
	}
	if (col)
	{
		std::reverse(temp.begin(), temp.end());
	}

	for (int m = 0; m < temp.size() - 1; ++m)
	{
		if (temp[m] == temp[m + 1])
		{
			temp[m] *= 2;
			temp[m + 1] = 0;
			m++;
		}
	}

	for (int m = 0; m < SIZE; ++m)
	{
		grid[row][m]->setValue(0);
		grid[row][m]->draw();
	}

	int tempj = col ? (SIZE-1) : 0;
	for (int m = 0; m < temp.size(); ++m)
	{
		if (temp[m] != 0)
		{
			grid[row][tempj]->setValue(temp[m]);
			grid[row][tempj]->draw();
			if (col)
				tempj--;
			else
				tempj++;
		}
	}
}

void Game::mergeFunc(bool row, int col)
{
	QList<int> temp;
	for (int i = 0; i != SIZE; ++i)
	{
		if (grid[i][col]->getValue() != 0)
		{
			temp.append(grid[i][col]->getValue());
		}
	}
	if (row)
	{
		std::reverse(temp.begin(), temp.end());
	}

	for (int m = 0; m < temp.size() - 1; ++m)
	{
		if (temp[m] == temp[m + 1])
		{
			temp[m] *= 2;
			temp[m + 1] = 0;
			m++;
		}
	}

	for (int m = 0; m < SIZE; ++m)
	{
		grid[m][col]->setValue(0);
		grid[m][col]->draw();
	}

	int tempj = row ? (SIZE-1) : 0;
	for (int m = 0; m < temp.size(); ++m)
	{
		if (temp[m] != 0)
		{
			grid[tempj][col]->setValue(temp[m]);
			grid[tempj][col]->draw();
			if (row)
				tempj--;
			else
				tempj++;
		}
	}
}

void Game::keyReleaseEvent(QKeyEvent* event)
{
	switch (event->key())
	{
	case Qt::Key_Up:
		this->merge(0);
		break;
	case Qt::Key_Down:
		this->merge(1);
		break;
	case Qt::Key_Left:
		this->merge(2);
		break;
	case Qt::Key_Right:
		this->merge(3);
		break;
	case Qt::Key_Escape:
		for (int i = 0; i != SIZE; ++i)
		{
			for (int j = 0; j != SIZE; ++j)
			{
				grid[i][j]->setValue(0);
				grid[i][j]->draw();
			}
		}
		break;
	default:
		break;
	}
}
