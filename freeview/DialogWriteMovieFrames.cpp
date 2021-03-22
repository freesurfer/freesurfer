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
#include "DialogWriteMovieFrames.h"
#include "ui_DialogWriteMovieFrames.h"
#include "MainWindow.h"
#include "RenderView.h"
#include "RenderView2D.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QCloseEvent>
#include "Layer.h"
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <QDebug>
#include <QSettings>

DialogWriteMovieFrames::DialogWriteMovieFrames(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogWriteMovieFrames)
{
  ui->setupUi(this);
  ui->doubleSpinBoxStep->hide();
  ui->labelSliceStepDouble->hide();
  connect(MainWindow::GetMainWindow(), SIGNAL(MainViewChanged(int)),
          this, SLOT(UpdateUI()));
  m_timer.setInterval(1000);
  connect(&m_timer, SIGNAL(timeout()), this, SLOT(OnTimeOut()));
  QSettings s;
  ui->lineEditOutputLocation->setText(s.value("MovieFrames/OutputLocation").toString());
}

DialogWriteMovieFrames::~DialogWriteMovieFrames()
{
  QSettings s;
  s.setValue("MovieFrames/OutputLocation", ui->lineEditOutputLocation->text().trimmed());
  delete ui;
}

void DialogWriteMovieFrames::showEvent(QShowEvent *e)
{
  UpdateUI();

  QDialog::showEvent(e);
}

void DialogWriteMovieFrames::closeEvent(QCloseEvent *e)
{
  if (m_timer.isActive())
  {
    if ( QMessageBox::question(this, "Close", "Still writing movie frames. Do you want to abort?") != QMessageBox::Yes)
    {
      e->ignore();
      return;
    }
  }

  QDialog::closeEvent(e);
}

void DialogWriteMovieFrames::UpdateUI(bool bUpdateNumbers)
{
  if (bUpdateNumbers)
  {
    m_view = MainWindow::GetMainWindow()->GetMainView();
    m_b3D = (MainWindow::GetMainWindow()->GetMainViewId() == 3);
    ui->comboBoxFlyThrough->setItemData(2, (m_b3D?33:0), Qt::UserRole - 1);
    if (!m_b3D && ui->comboBoxFlyThrough->currentIndex() == 2)
      ui->comboBoxFlyThrough->setCurrentIndex(0);
    LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer("MRI");
    if ( mri )
    {
      ui->comboBoxFlyThrough->setItemData(1, (mri->GetNumberOfFrames()>1?33:0), Qt::UserRole-1);
    }
  }
  QList<QWidget*> widgets = findChildren<QWidget*>();
  foreach (QWidget* w, widgets)
  {
    if (w != ui->pushButtonAbort && w != ui->pushButtonClose && w != ui->pushButtonWrite)
      w->setEnabled(!m_timer.isActive());
  }

  ui->pushButtonAbort->setEnabled(m_timer.isActive());
  ui->pushButtonClose->setEnabled(!m_timer.isActive());
  ui->pushButtonWrite->setEnabled(!m_timer.isActive());
}

void DialogWriteMovieFrames::OnComboBoxFlyThrough(int nIndex)
{
  int nPlane = MainWindow::GetMainWindow()->GetMainViewId();
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer("MRI");
  if ( mri )
  {
    if (nIndex == 0)
    {
      vtkImageData* imagedata = mri->GetImageData();
      int* dim = imagedata->GetDimensions();
      ui->spinBoxEnd->setValue(dim[nPlane]-1);
    }
    else if (nIndex == 1)
    {
      ui->spinBoxEnd->setValue(mri->GetNumberOfFrames()-1);
    }
  }
  if (nIndex == 2)
  {
    ui->spinBoxEnd->setValue(360);
  }
  ui->spinBoxStart->setValue(0);
  ui->doubleSpinBoxStep->setVisible(nIndex == 2);
  ui->labelSliceStepDouble->setVisible(nIndex == 2);
  ui->labelSliceStep->setVisible(nIndex != 2);
  ui->spinBoxStep->setVisible(nIndex != 2);
}

void DialogWriteMovieFrames::OnOpen()
{
  QString strg = QFileDialog::getExistingDirectory(this, "Select Output Directory");
  if (!strg.isEmpty())
  {
    ui->lineEditOutputLocation->setText(strg);
    ui->lineEditOutputLocation->setCursorPosition(ui->lineEditOutputLocation->text().length());
  }
}

void DialogWriteMovieFrames::OnWrite()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  SettingsScreenshot ss = mainwnd->GetScreenShotSettings();
  ss.HideCoords = !ui->checkBoxShowAnnotations->isChecked();
  mainwnd->SetScreenShotSettings(ss);

  m_strOutputDir = ui->lineEditOutputLocation->text().trimmed();
  if (m_strOutputDir.isEmpty())
  {
    QMessageBox::warning(this, "Error", "Please enter output directory to save movie frames.");
    return;
  }
  QDir dir(m_strOutputDir);
  if (!dir.exists())
  {
    QMessageBox::warning(this, "Error", "Output directory does not exist.");
    return;
  }
  if (ui->spinBoxStep->value() == 0)
  {
    QMessageBox::warning(this, "Error", "Step size can not be 0.");
    return;
  }
  m_strPrefix = ui->lineEditFilenamePrefix->text().trimmed();
  if (m_strOutputDir[m_strOutputDir.length()-1] != '/' &&
      m_strOutputDir[m_strOutputDir.length()-1] != '\\')
  {
    m_strOutputDir += QDir::separator();
  }
  m_nStepSize = ui->spinBoxStep->value();
  m_dStepSize = ui->doubleSpinBoxStep->value();
  m_nStartNumber = ui->spinBoxStart->value();
  m_nStepCount = 0;
  m_nTotalSteps = 1;
  int nIndex = ui->comboBoxFlyThrough->currentIndex();
  if (nIndex == 0)    // slice
  {
    MainWindow* mwnd = MainWindow::GetMainWindow();
    Layer* layer = mwnd->GetActiveLayer( "MRI" );
    if (!layer)
    {
      layer = mwnd->GetActiveLayer("Surface");
    }
    if (layer)
    {
      m_nTotalSteps = (ui->spinBoxEnd->value()-ui->spinBoxStart->value())/m_nStepSize+1;
    }
    ((RenderView2D*)m_view)->SetSliceNumber( m_nStartNumber );
  }
  else if (nIndex == 1) // frame
  {
    MainWindow* mwnd = MainWindow::GetMainWindow();
    Layer* layer = mwnd->GetActiveLayer( "MRI" );
    if (layer)
    {
      m_nTotalSteps = (ui->spinBoxEnd->value()-ui->spinBoxStart->value())/m_nStepSize+1;
      ((LayerMRI*)layer)->SetActiveFrame(m_nStartNumber);
    }
  }
  else // angle
  {
    m_nTotalSteps = (int)((ui->spinBoxEnd->value()-ui->spinBoxStart->value())/m_dStepSize+1);
  }

  if (m_nTotalSteps < 1)
  {
    m_nTotalSteps = 1;
  }
  m_timer.start();
  UpdateUI(false);
  emit Started();
}

void DialogWriteMovieFrames::OnAbort()
{
  m_timer.stop();
  UpdateUI(false);
  emit Stopped();
}

void DialogWriteMovieFrames::OnTimeOut()
{
  QString fn;
  SettingsScreenshot settings = MainWindow::GetMainWindow()->GetScreenShotSettings();
  int nIndex = ui->comboBoxFlyThrough->currentIndex();
  int nFieldWidth = qMax(3, QString::number(m_nTotalSteps).size());
  if (nIndex == 0)    // slice
  {
    int nStart = m_nStartNumber+m_nStepSize*m_nStepCount;
    fn = QString("%1%2%3.%4").arg(m_strOutputDir)
        .arg(m_strPrefix)
        .arg(nStart, nFieldWidth, 10, QLatin1Char('0'))
        .arg(ui->comboBoxExtension->currentText());
    m_view->SaveScreenShot( fn,settings.AntiAliasing, settings.Magnification );
    if ( m_nStepCount+1 >= m_nTotalSteps || !((RenderView2D*)m_view)->SetSliceNumber( nStart + m_nStepSize ) )
    {
      OnAbort();
    }
  }
  else if (nIndex == 1)  // frame
  {
    int nStart = m_nStartNumber+m_nStepSize*m_nStepCount;
    fn = QString("%1%2%3.%4").arg(m_strOutputDir)
        .arg(m_strPrefix)
        .arg(nStart, nFieldWidth, 10, QLatin1Char('0'))
        .arg(ui->comboBoxExtension->currentText());
    m_view->SaveScreenShot( fn,settings.AntiAliasing, settings.Magnification );
    LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer( "MRI" ));
    if (m_nStepCount+1 < m_nTotalSteps && layer)
      layer->SetActiveFrame(nStart + m_nStepSize);
  }
  else          // angle
  {
    fn = QString("%1%2%3.%4").arg(m_strOutputDir)
        .arg(m_strPrefix)
        .arg(m_nStepCount, nFieldWidth, 10, QLatin1Char('0'))
        .arg(ui->comboBoxExtension->currentText());
    m_view->SaveScreenShot( fn, settings.AntiAliasing, settings.Magnification );
    CameraOperations ops;
    ops << CameraOperation("azimuth", m_dStepSize);
    m_view->SetCameraOperations(ops);
  }

  m_view->RequestRedraw(true);
  m_nStepCount++;
  if (m_nStepCount >= m_nTotalSteps)
  {
    OnAbort();
  }
  else
  {
    emit Progress(m_nStepCount*100/m_nTotalSteps);
  }
}
