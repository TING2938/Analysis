#include "Square.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	Square w;
	w.show();
	return a.exec();
}
