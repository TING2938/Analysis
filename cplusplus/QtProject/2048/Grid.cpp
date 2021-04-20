#include "Grid.h"
#include <QGraphicsDropShadowEffect>

Grid::Grid(QWidget *parent)
	: QLabel(parent)
{
	this->setFixedSize(110, 110);
    this->setAlignment(Qt::AlignCenter);
}

Grid::~Grid()
{
}

void Grid::setValue(int value)
{
	this->num = value;
}

int Grid::getValue() const
{
	return num;
}

void Grid::draw()
{
    if (num == 0)
    {
        setText("");
        setStyleSheet("Grid { background: rgb(204,192,179); border-radius: 10px; }");
    }
    else
    {
        setText(QString::number(num));
        switch (num)
        {
        case 2: {
            setStyleSheet("Grid { background: rgb(238,228,218); color: rgb(119,110,101); font: bold; border-radius: 10px; font: 40pt; }");
            break;
        }
        case 4: {
            setStyleSheet("Grid { background: rgb(237,224,200); color: rgb(119,110,101); font: bold; border-radius: 10px; font: 40pt; }");
            break;
        }
        case 8: {
            setStyleSheet("Grid { background: rgb(242,177,121); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 40pt; }");
            break;
        }
        case 16: {
            setStyleSheet("Grid { background: rgb(245,150,100); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 40pt; }");
            break;
        }
        case 32: {
            setStyleSheet("Grid { background: rgb(245,125,95); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 40pt; }");
            break;
        }
        case 64: {
            setStyleSheet("Grid { background: rgb(245,95,60); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 40pt; }");
            break;
        }
        case 128: {
            setStyleSheet("Grid { background: rgb(237,207,114); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 30pt; }");
            break;
        }
        case 256: {
            setStyleSheet("Grid { background: rgb(237,204,97); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 30pt; }");
            break;
        }
        case 512: {
            setStyleSheet("Grid { background: rgb(237,204,97); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 30pt; }");
            break;
        }
        case 1024: {
            setStyleSheet("Grid { background: rgb(237,204,97); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 30pt; }");
            break;
        }
        case 2048: {
            setStyleSheet("Grid { background: rgb(237,204,97); color: rgb(255,255,255); font: bold; border-radius: 10px; font: 20pt; }");
            break;
        }
        default: {
            setStyleSheet("Grid { background: rgb(238,228,218); color: rgb(119,110,101); font: bold; border-radius: 10px; font: 20pt; }");
            break;
        }
        }
    }
}
