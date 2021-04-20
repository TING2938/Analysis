#include "Grid.h"
#include <random>
#include <QGridLayout>

Grid::Grid(QWidget *parent)
	: QWidget(parent)
{
	a = 25;
	margin = QPoint(3*a, 3*a);
	nrow = 20;
	ncol = 10;
	ntype = 8;

	palette.resize(ntype);
	for (int i = 0; i != ntype; ++i)
	{
		palette[i].setColor(QPalette::Background, QColor(255*i/ntype, 150, 150));
	}

	map.resize(nrow);
	for (auto&& m : map)
	{
		m.resize(ncol);
		m.fill(0);
	}

	labelMap.resize(nrow);
	for (int i = 0; i != nrow; ++i)
	{
		labelMap[i].resize(ncol);
		for (int j = 0; j != ncol; ++j)
		{
			labelMap[i][j] = new QLabel(parent);
			labelMap[i][j]->resize(a, a);
			labelMap[i][j]->move(margin + a * QPoint(j, i));
			labelMap[i][j]->setFrameShape(QFrame::Box);
			labelMap[i][j]->setAutoFillBackground(true);
			labelMap[i][j]->hide();
		}
	}

	setType(getRand(0, ntype-1));

	for (int i = 0; i < 4; ++i)
	{
		grid[i] = new QLabel(parent);
		grid[i]->resize(a, a);
		grid[i]->setFrameShape(QFrame::Box); 
		grid[i]->setAutoFillBackground(true);
		grid[i]->setPalette(palette[type]);

		pos[i] = typeVec[i] + QPoint(ncol/2-1, 0);
	}
	setPos();
}

Grid::~Grid()
{
}


void Grid::setType(int t)
{
	type = t;
	switch (type)
	{
		// ����
	case 0:
		typeVec[0] = QPoint(0, 0);
		typeVec[1] = QPoint(1, 0);
		typeVec[2] = QPoint(2, 0);
		typeVec[3] = QPoint(3, 0);
		break;

		// T��
	case 1:
		typeVec[0] = QPoint(0, 0);
		typeVec[1] = QPoint(1, 0);
		typeVec[2] = QPoint(2, 0);
		typeVec[3] = QPoint(1, 1);
		break;

		// ����
	case 2:
		typeVec[0] = QPoint(0, 0);
		typeVec[1] = QPoint(0, 1);
		typeVec[2] = QPoint(0, 2);
		typeVec[3] = QPoint(0, 3);
		break;

		// ����
	case 3:
		typeVec[0] = QPoint(0, 0);
		typeVec[1] = QPoint(0, 1);
		typeVec[2] = QPoint(1, 0);
		typeVec[3] = QPoint(1, 1);
		break;

		// S��
	case 4:
		typeVec[0] = QPoint(0, 0);
		typeVec[1] = QPoint(0, 1);
		typeVec[2] = QPoint(1, 1);
		typeVec[3] = QPoint(1, 2);
		break;

		// L��
	case 5:
		typeVec[0] = QPoint(0, 0);
		typeVec[1] = QPoint(0, 1);
		typeVec[2] = QPoint(0, 2);
		typeVec[3] = QPoint(1, 2);
		break;

		// ��L��
	case 6:
		typeVec[0] = QPoint(1, 0);
		typeVec[1] = QPoint(1, 1);
		typeVec[2] = QPoint(1, 2);
		typeVec[3] = QPoint(0, 2);
		break;

		// Z��
	case 7:
		typeVec[0] = QPoint(0, 0);
		typeVec[1] = QPoint(1, 0);
		typeVec[2] = QPoint(1, 1);
		typeVec[3] = QPoint(2, 1);
		break;

	default:
		break;
	}

}

void Grid::rotate(int dir)
{
	int SIN = dir;

	if (type != 3)
	{
		QPoint tmp[4];
		for (int i = 0; i != 4; ++i)
		{
			tmp[i].rx() = - SIN * (pos[i].y() - pos[1].y()) + pos[1].x();
			tmp[i].ry() = SIN * (pos[i].x() - pos[1].x()) + pos[1].y();

			if (tmp[i].x() < 0 || tmp[i].x() >= ncol ||
				tmp[i].y() < 0 || tmp[i].y() >= nrow ||
				map[tmp[i].y()][tmp[i].x()] == 1)
			{
				return;
			}
		}
		for (int i = 0; i != 4; ++i)
		{
			pos[i] = tmp[i];
		}
		setPos();
	}
}

void Grid::moveLeft()
{
	for (auto&& p : pos)
	{
		if (p.x() == 0 || map[p.y()][p.x() -1] == 1)
			return;
	}

	for (auto&& p : pos)
	{
		p.rx() -= 1;
	}
	setPos();
}

void Grid::moveRight()
{
	for (auto&& p : pos)
	{
		if (p.x() == ncol-1 || map[p.y()][p.x()+1] == 1)
			return;
	}

	for (auto&& p : pos)
	{
		p.rx() += 1;
	}
	setPos();
}

void Grid::drop()
{
	for (int i = 0; i < 4; ++i)
	{
		if (pos[i].y() == nrow-1 || map[pos[i].y()+1][pos[i].x()] == 1)
		{
			for (auto&& p : pos)
			{
				map[p.y()][p.x()] = 1;
				labelMap[p.y()][p.x()]->setPalette(palette[type]);
			}
			checkMap();
			drawMap();
			this->reset();
			return;
		}
	}
	
	for (auto&& p : pos)
	{
		p.ry() += 1;
	}
	setPos();
}

void Grid::reset()
{
	setType(getRand(0, ntype-1));

	for (int i = 0; i < 4; ++i)
	{
		pos[i] = typeVec[i] + QPoint(ncol/2-1, 0);
		grid[i]->setPalette(palette[type]);
	}
	setPos();
}

void Grid::drawMap()
{
	for (int i = 0; i < nrow; ++i)
	{
		for (int j = 0; j < ncol; ++j)
		{
			if (map[i][j])
			{
				labelMap[i][j]->show();
			}
			else
			{
				labelMap[i][j]->hide();
			}
		}
	}
}

void Grid::checkMap()
{
	for (int i = 0; i != nrow; ++i)
	{
		if (!map[i].contains(0))
		{
			for (int k = i; k >= 1; --k)
			{
				map[k] = map[k - 1];
			}
		}
	}
}

void Grid::setPos()
{
	for (int i = 0; i < 4; ++i)
	{
		grid[i]->move(margin + a * pos[i]);
	}
}

int Grid::getRand(int beg, int end)
{
	static std::random_device rd;
	static std::default_random_engine e(rd());
	std::uniform_int_distribution<int> dist(beg, end);
	return dist(e);
}

