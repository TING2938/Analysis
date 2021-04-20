#pragma once

#include <QtWidgets/QMainWindow>
#include <QTextEdit>


class Notepad : public QMainWindow
{
	Q_OBJECT

public:
	Notepad(QWidget *parent = Q_NULLPTR);


	void saveFile();

private:
	QMenu* menu_file;
	QMenu* menu_edit;
	QAction* act_new;
	QAction* act_save;
	QAction* act_open;
	QAction* act_edit;
	QTextEdit* text;
};
