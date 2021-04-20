#pragma once

#include <QLabel>

class Grid : public QLabel
{
	Q_OBJECT

public:
	Grid(QWidget *parent= nullptr);
	~Grid();

	void setValue(int value);
	int getValue() const;

	void draw();
private:
	int num;
};
