/**
 * @file  WindowTimeCourse.cpp
 * @brief Tool window to display time course data
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/08/29 15:25:00 $
 *    $Revision: 1.2 $
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

#include "WindowTimeCourse.h"
#include "ui_WindowTimeCourse.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerCollection.h"
#include "FSVolume.h"
#include <QSettings>

WindowTimeCourse::WindowTimeCourse(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WindowTimeCourse)
{
    ui->setupUi(this);
    this->setWindowFlags(Qt::Tool);
    this->setWindowTitle("Time Course");
    connect(ui->widgetPlot, SIGNAL(FrameChanged(int)), this, SLOT(OnFrameChanged(int)));

    QSettings s;
    QVariant v = s.value("WindowTimeCourse/Geomerty");
    if (v.isValid())
      this->restoreGeometry(v.toByteArray());
    ui->checkBoxAutoScale->setChecked(s.value("WindowTimeCourse/AutoScale", true).toBool());
}

WindowTimeCourse::~WindowTimeCourse()
{
  QSettings s;
  s.setValue("WindowTimeCourse/Geomerty", this->saveGeometry());
  s.setValue("WindowTimeCourse/AutoScale", ui->checkBoxAutoScale->isChecked());
  delete ui;
}

void WindowTimeCourse::UpdateData()
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer("MRI"));
  if (layer && layer->GetNumberOfFrames() > 1)
  {
    double ras[3];
    int n[3];
    MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetSlicePosition(ras);
    layer->RemapPositionToRealRAS(ras, ras);
    layer->RASToOriginalIndex(ras, n);
    QList<double> data;
    for (int i = 0; i < layer->GetNumberOfFrames(); i++)
      data <<  layer->GetVoxelValueByOriginalIndex(n[0], n[1], n[2], i);
    FSVolume* vol = layer->GetSourceVolume();
    ui->widgetPlot->SetTimeCourseData(data, vol->GetMinValue(), vol->GetMaxValue(), layer->GetTR());
    ui->widgetPlot->SetCurrentFrame(layer->GetActiveFrame());
  }
}

void WindowTimeCourse::OnFrameChanged(int frame)
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer("MRI"));
  if (layer && frame != layer->GetActiveFrame() && frame < layer->GetNumberOfFrames())
  {
    layer->SetActiveFrame(frame);
  }
}

void WindowTimeCourse::SetCurrentFrame(int n)
{
  ui->widgetPlot->SetCurrentFrame(n);
}
