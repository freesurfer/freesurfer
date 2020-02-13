#include "DialogNewAnnotation.h"
#include "ui_DialogNewAnnotation.h"
#include <QFileDialog>
#include <QMessageBox>

DialogNewAnnotation::DialogNewAnnotation(QWidget *parent, const QString& dir) :
  QDialog(parent), m_strDir(dir),
  ui(new Ui::DialogNewAnnotation)
{
  ui->setupUi(this);
}

DialogNewAnnotation::~DialogNewAnnotation()
{
  delete ui;
}

void DialogNewAnnotation::OnOK()
{
  if (GetName().isEmpty())
  {
    QMessageBox::information(this, "Error", "Please enter a name for the new annotation");
    return;
  }
//  else if (GetColorTableFile().isEmpty())
//  {
//    QMessageBox::information(this, "Error", "Please select a color table file for the new annotation");
//    return;
//  }
  accept();
}

void DialogNewAnnotation::OnOpen()
{
  QString fn = QFileDialog::getOpenFileName( this, "Select Color Table file",
                                     m_strDir,
                                     "Color Table files (*.txt *.ctab)");
  if (!fn.isEmpty())
  {
    ui->lineEditColorTableFile->setText(fn);
  }
}

QString DialogNewAnnotation::GetName()
{
  return ui->lineEditName->text().trimmed();
}

QString DialogNewAnnotation::GetColorTableFile()
{
  return ui->lineEditColorTableFile->text().trimmed();
}
