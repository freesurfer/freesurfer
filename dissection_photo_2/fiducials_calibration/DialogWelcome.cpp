#include "DialogWelcome.h"
#include "ui_DialogWelcome.h"

DialogWelcome::DialogWelcome(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogWelcome)
{
  ui->setupUi(this);
}

DialogWelcome::~DialogWelcome()
{
  delete ui;
}
