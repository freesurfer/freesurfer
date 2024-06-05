#include "ProgressWindow.h"
#include "ui_ProgressWindow.h"

ProgressWindow::ProgressWindow(QWidget *parent) :
  QFrame(parent),
  ui(new Ui::ProgressWindow)
{
  ui->setupUi(this);
  setWindowFlags(Qt::Tool | Qt::CustomizeWindowHint);
}

ProgressWindow::~ProgressWindow()
{
  delete ui;
}
