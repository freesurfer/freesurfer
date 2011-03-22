/**
 * @file  DialogSaveScreenshot.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/22 21:21:26 $
 *    $Revision: 1.9 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "DialogSaveScreenshot.h"
#include "ui_DialogSaveScreenshot.h"
#include "MainWindow.h"
#include "RenderView.h"
#include "MyUtils.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QFileInfo>
#include <QSettings>

DialogSaveScreenshot::DialogSaveScreenshot(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogSaveScreenshot)
{
  ui->setupUi(this);
  QSettings settings;
  ui->lineEditFileName->setText(settings.value("ScreenShot/LastSavedFile").toString());
}

DialogSaveScreenshot::~DialogSaveScreenshot()
{
  QSettings settings;
  settings.setValue("ScreenShot/LastSavedFile", GetFileName());
  delete ui;
}


QString DialogSaveScreenshot::GetFileName()
{
  QString filename = MyUtils::CygwinPathProof(ui->lineEditFileName->text().trimmed());;
  return QFileInfo(QDir::current(), filename).absoluteFilePath();
}

void DialogSaveScreenshot::SetSettings( SettingsScreenshot s )
{
  ui->checkBoxAntiAliasing->setChecked( s.AntiAliasing );
  ui->checkBoxHideCursor->setChecked( s.HideCursor );
  ui->checkBoxHideAnnotation->setChecked( s.HideCoords );
  ui->spinBoxMagnification->setValue( s.Magnification );
}

SettingsScreenshot DialogSaveScreenshot::GetSettings()
{
  SettingsScreenshot s;
  s.AntiAliasing  = ui->checkBoxAntiAliasing->isChecked();
  s.HideCursor    = ui->checkBoxHideCursor->isChecked();
  s.HideCoords    = ui->checkBoxHideAnnotation->isChecked();
  s.Magnification = ui->spinBoxMagnification->value();

  return s;
}

void DialogSaveScreenshot::OnOpen()
{
  QString dir = m_strLastDir;
  if (!GetFileName().isEmpty())
  {
    dir = QFileInfo(GetFileName()).absolutePath();
  }
  QString fn = QFileDialog::getSaveFileName( this, "Save Screenshot", dir, "All Files (*.*)" );
  if ( !fn.isEmpty() )
  {
    ui->lineEditFileName->setText(MyUtils::Win32PathProof(fn));
    ui->lineEditFileName->setCursorPosition(ui->lineEditFileName->text().size());
  }
}

void DialogSaveScreenshot::OnSave()
{
  if ( GetFileName().isEmpty() )
  {
    QMessageBox::warning(this, "Error", "Please enter file name to be saved.");
    return;
  }

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  mainwnd->SetScreenShotSettings(GetSettings());
  if (!mainwnd->GetMainView()->
      SaveScreenShot(GetFileName(), ui->checkBoxAntiAliasing, ui->spinBoxMagnification->value()))
  {
    QMessageBox::warning(this, "Error", "Failed to save screenshot. Please make sure the directory exists and writable.");
    return;
  }

  if (!ui->checkBoxKeepWindow->isChecked())
  {
    hide();
  }
}
