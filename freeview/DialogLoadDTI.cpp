/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "DialogLoadDTI.h"
#include "ui_DialogLoadDTI.h"
#include "MainWindow.h"
#include <QMessageBox>
#include <QFileInfo>
#include <QFileDialog>
#include "MyUtils.h"

DialogLoadDTI::DialogLoadDTI(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLoadDTI)
{
  ui->setupUi(this);
}

DialogLoadDTI::~DialogLoadDTI()
{
  delete ui;
}

QString DialogLoadDTI::GetVectorFileName()
{
  return MyUtils::CygwinPathProof(ui->lineEditVector->text().trimmed());
}

QString DialogLoadDTI::GetFAFileName()
{
  return MyUtils::CygwinPathProof(ui->lineEditFA->text().trimmed());
}

QString DialogLoadDTI::GetEigenvalueFileName()
{
  return MyUtils::CygwinPathProof(ui->lineEditEigenvalue->text().trimmed());
}

QString DialogLoadDTI::GetRegFileName()
{
  if ( ui->checkBoxRegistration->isChecked() )
  {
    return MyUtils::CygwinPathProof(ui->lineEditRegistration->text().trimmed());
  }
  else
  {
    return "";
  }
}

void DialogLoadDTI::OnOK()
{
  if ( GetVectorFileName().isEmpty() )
  {
    QMessageBox::warning(
          this, "Error",
          "Vector file name can not be empty.");
    return;
  }
  else if ( GetFAFileName().isEmpty() )
  {
    QMessageBox::warning(
          this, "Error",
          "FA file name can not be empty.");
    return;
  }
  else if ( ui->checkBoxRegistration->isChecked() && GetRegFileName().isEmpty() )
  {
    QMessageBox::warning(
          this, "Error",
          "Registration file name can not be empty.");
    return;
  }

  accept();
}

void DialogLoadDTI::OnButtonVector()
{
  QString filename = QFileDialog::getOpenFileName(
        this,
        "Select vector file",
        MainWindow::AutoSelectLastDir( m_strLastDir, "mri" ),
        "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditVector->setText( MyUtils::Win32PathProof(filename) );
    ui->lineEditVector->setCursorPosition( ui->lineEditVector->text().size() );
    m_strLastDir = QFileInfo(filename).canonicalPath();
  }
}

void DialogLoadDTI::OnButtonFA()
{
  QString filename = QFileDialog::getOpenFileName(
        this,
        "Select FA file",
        MainWindow::AutoSelectLastDir( m_strLastDir, "mri" ),
        "FA files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditFA->setText( MyUtils::Win32PathProof(filename) );
    ui->lineEditFA->setCursorPosition( ui->lineEditFA->text().size() );
    m_strLastDir = QFileInfo(filename).canonicalPath();
  }
}

void DialogLoadDTI::OnButtonRegistration()
{
  QString filename = QFileDialog::getOpenFileName(
        this,
        "Select registration file",
        MainWindow::AutoSelectLastDir( m_strLastDir, "mri" ),
        "Registration files (*)" );
  if ( !filename.isEmpty() )
  {
    ui->lineEditRegistration->setText( MyUtils::Win32PathProof(filename) );
    ui->lineEditRegistration->setCursorPosition( ui->lineEditRegistration->text().size() );
  }
}

void DialogLoadDTI::OnButtonEigenvalue()
{
  QString filename = QFileDialog::getOpenFileName(
        this,
        "Select eigenvalue file",
        MainWindow::AutoSelectLastDir( m_strLastDir, "mri" ),
        "FA files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditEigenvalue->setText( MyUtils::Win32PathProof(filename) );
    ui->lineEditEigenvalue->setCursorPosition( ui->lineEditEigenvalue->text().size() );
    m_strLastDir = QFileInfo(filename).canonicalPath();
  }
}


bool DialogLoadDTI::IsToResample()
{
  return ui->checkBoxResample->isChecked();
}


void DialogLoadDTI::Initialize( bool bResample, bool bEnableCheckBox )
{
  ui->checkBoxResample->setChecked( bResample );
  ui->checkBoxResample->setEnabled( bEnableCheckBox );
}
