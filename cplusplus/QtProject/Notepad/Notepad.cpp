#include "Notepad.h"
#include <QMenuBar>
#include <QMenu>
#include <QTextEdit>
#include <QFile>
#include <QString>
#include <QFileDialog>

Notepad::Notepad(QWidget *parent)
	: QMainWindow(parent)
{
	this->resize(1600, 1000);

	auto bar = menuBar();
	menu_file = bar->addMenu("File");
	menu_edit = bar->addMenu("Edit");
	act_new = menu_file->addAction("New");
	act_open = menu_file->addAction("Open");
	act_save = menu_file->addAction("Save");
	auto act_exit = menu_file->addAction("Exit");
	act_edit = menu_edit->addAction("Edit");
	
	text = new QTextEdit(this);
	this->setCentralWidget(text);
	text->setFrameShape(QFrame::NoFrame);
	
	connect(act_exit, &QAction::triggered, this, &QWidget::close);
	connect(act_save, &QAction::triggered, this, &Notepad::saveFile);
}

void Notepad::saveFile()
{
	QString str = text->toPlainText();
	auto dir = QDir::toNativeSeparators(QFileDialog::getSaveFileName(this, "Save path"));
	QFile file(dir);
	file.open(QIODevice::WriteOnly);
	file.write(str.toUtf8());
	file.close();
}
