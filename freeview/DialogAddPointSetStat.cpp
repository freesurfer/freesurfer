#include "DialogAddPointSetStat.h"
#include "ui_DialogAddPointSetStat.h"
#include <QMessageBox>

DialogAddPointSetStat::DialogAddPointSetStat(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogAddPointSetStat)
{
  ui->setupUi(this);
}

DialogAddPointSetStat::~DialogAddPointSetStat()
{
  delete ui;
}

void DialogAddPointSetStat::OnButtonAdd()
{
  bool bOK;
  ui->lineEditValue->text().toDouble(&bOK);
  if (ui->lineEditName->text().trimmed().isEmpty())
  {
    QMessageBox::warning(this, "Error", "Please enter a name");
  }
  else if (!bOK)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid value");
  }
  else
    accept();
}

QString DialogAddPointSetStat::GetStatName()
{
  return ui->lineEditName->text().trimmed();
}

double DialogAddPointSetStat::GetStatValue()
{
  return ui->lineEditValue->text().toDouble();
}
