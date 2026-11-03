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

void ProgressWindow::SetIndeterminate(const QString& msg)
{
  ui->label->setText(msg);
  ui->progressBar->setTextVisible(false);
  ui->progressBar->setMinimum(0);
  ui->progressBar->setMaximum(0);
  ui->progressBar->setValue(0);
}
