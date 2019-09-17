/**
 * @file  DialogSaveScreenshot.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2015/09/16 20:36:43 $
 *    $Revision: 1.10 $
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
#include <QDebug>

DialogSaveScreenshot::DialogSaveScreenshot(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogSaveScreenshot)
{
  ui->setupUi(this);
  //  QSettings settings;
  //  ui->lineEditFileName->setText(settings.value("ScreenShot/LastSavedFile").toString());
  m_strLastDir = QDir::currentPath();
#ifdef Q_OS_MAC
  ui->checkBoxAutoTrim->hide();
#endif
//  qDebug() << m_strLastDir;
}

DialogSaveScreenshot::~DialogSaveScreenshot()
{
  //  QSettings settings;
  //  settings.setValue("ScreenShot/LastSavedFile", GetFileName());
  delete ui;
}


QString DialogSaveScreenshot::GetFileName()
{
  QString filename = ui->lineEditFileName->text().trimmed();
  if (!filename.isEmpty())
    return QFileInfo(QDir::current(), filename).absoluteFilePath();
  else
    return "";
}

void DialogSaveScreenshot::SetSettings( SettingsScreenshot s )
{
  ui->checkBoxAntiAliasing->setChecked( s.AntiAliasing );
  ui->checkBoxHideCursor->setChecked( s.HideCursor );
  ui->checkBoxHideAnnotation->setChecked( s.HideCoords );
  ui->spinBoxMagnification->setValue( s.Magnification );
  ui->checkBoxAutoTrim->setChecked( s.AutoTrim );
  ui->checkBoxHideScaleBar->setChecked( s.HideScaleBar );
}

SettingsScreenshot DialogSaveScreenshot::GetSettings()
{
  SettingsScreenshot s;
  s.AntiAliasing  = ui->checkBoxAntiAliasing->isChecked();
  s.HideCursor    = ui->checkBoxHideCursor->isChecked();
  s.HideCoords    = ui->checkBoxHideAnnotation->isChecked();
  s.Magnification = ui->spinBoxMagnification->value();
  s.AutoTrim  = ui->checkBoxAutoTrim->isChecked();
  s.HideScaleBar  = ui->checkBoxHideScaleBar->isChecked();

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
      SaveScreenShot(GetFileName(), ui->checkBoxAntiAliasing, ui->spinBoxMagnification->value(), ui->checkBoxAutoTrim->isChecked()))
  {
    QMessageBox::warning(this, "Error", "Failed to save screenshot. Please make sure the directory exists and writable.");
    return;
  }
  m_strLastDir = QFileInfo(GetFileName()).absolutePath();

  if (!ui->checkBoxKeepWindow->isChecked())
  {
    hide();
  }
}
