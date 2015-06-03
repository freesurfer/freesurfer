#include "DialogThresholdFilter.h"
#include "ui_DialogThresholdFilter.h"

DialogThresholdFilter::DialogThresholdFilter(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogThresholdFilter)
{
    ui->setupUi(this);
}

DialogThresholdFilter::~DialogThresholdFilter()
{
    delete ui;
}

double DialogThresholdFilter::GetInValue()
{
  return ui->lineEditInValue->text().toDouble();
}

double DialogThresholdFilter::GetOutValue()
{
  return ui->lineEditOutValue->text().toDouble();
}

bool DialogThresholdFilter::GetReplaceIn()
{
  return ui->checkBoxReplaceIn->isChecked();
}

bool DialogThresholdFilter::GetReplaceOut()
{
  return ui->checkBoxReplaceOut->isChecked();
}

void DialogThresholdFilter::GetThreshold(double *th)
{
  th[0] = ui->lineEditThresholdLow->text().toDouble();
  if (ui->lineEditThresholdHigh->text().toLower() == "inf")
    th[1] = 1e20;
  else
    th[1] = ui->lineEditThresholdHigh->text().toDouble();
}
