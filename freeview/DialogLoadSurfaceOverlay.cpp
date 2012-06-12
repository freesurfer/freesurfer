#include "DialogLoadSurfaceOverlay.h"
#include "ui_DialogLoadSurfaceOverlay.h"
#include <QMessageBox>
#include <QFileDialog>

DialogLoadSurfaceOverlay::DialogLoadSurfaceOverlay(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogLoadSurfaceOverlay)
{
    ui->setupUi(this);
}

DialogLoadSurfaceOverlay::~DialogLoadSurfaceOverlay()
{
    delete ui;
}

void DialogLoadSurfaceOverlay::OnOK()
{
  if (ui->checkBoxRegistration->isChecked() && GetRegistration().isEmpty())
  {
    QMessageBox::warning(this, "Error", "Please enter a valid registration file name");
  }
  else
    accept();
}

QString DialogLoadSurfaceOverlay::GetFileName()
{
  return ui->lineEditFile->text().trimmed();
}

QString DialogLoadSurfaceOverlay::GetRegistration()
{
  if (ui->checkBoxRegistration->isChecked())
    return ui->lineEditRegistration->text().trimmed();
  else
    return QString();
}

void DialogLoadSurfaceOverlay::OnButtonOpen()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select volume files",
                          m_strLastDir,
                          "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditFile->setText(filename);
    ui->lineEditFile->setCursorPosition( filename.size() );
  }
}

void DialogLoadSurfaceOverlay::OnButtonRegistration()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select registration file",
                     m_strLastDir,
                     "Registration files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditRegistration->setText(filename);
    ui->lineEditRegistration->setCursorPosition( filename.size() );
  }
}
