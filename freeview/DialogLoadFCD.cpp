#include "DialogLoadFCD.h"
#include "ui_DialogLoadFCD.h"

DialogLoadFCD::DialogLoadFCD(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLoadFCD)
{
  ui->setupUi(this);
}

DialogLoadFCD::~DialogLoadFCD()
{
  delete ui;
}
