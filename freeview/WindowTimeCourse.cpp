/**
 * @file  WindowTimeCourse.cpp
 * @brief Tool window to display time course data
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2014/02/27 21:05:45 $
 *    $Revision: 1.4 $
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
#include "LayerSurface.h"
#include "SurfaceOverlay.h"
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
  QString type = MainWindow::GetMainWindow()->GetCurrentLayerType();
  if (type == "MRI")
  {
    LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer(type));
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
      connect(layer, SIGNAL(CorrelationSurfaceChanged(LayerSurface*)),
              this, SLOT(OnLayerCorrelationSurfaceChanged()), Qt::UniqueConnection);
      setWindowTitle(QString("Time Course (%1)").arg(layer->GetName()));
    }
  }
  else if (type == "Surface")
  {
    LayerSurface* surf = qobject_cast<LayerSurface*>(MainWindow::GetMainWindow()->GetActiveLayer(type));
    if (surf && surf->GetActiveOverlay() && surf->GetActiveOverlay()->GetNumberOfFrames() > 1)
    {
      double pos[3];
      MainWindow::GetMainWindow()->GetLayerCollection("Surface")->GetSlicePosition(pos);
      int nVert = surf->GetVertexIndexAtTarget(pos, NULL);
      if (nVert < 0)
        return;

      SurfaceOverlay* overlay = surf->GetActiveOverlay();
      int nFrames = overlay->GetNumberOfFrames();
      float* buffer = new float[nFrames];
      double range[2];
      overlay->GetDataAtVertex(nVert, buffer);
      overlay->GetRawRange(range);
      QList<double> data;
      for (int i = 0; i < nFrames; i++)
        data << buffer[i];
      delete[] buffer;
      ui->widgetPlot->SetTimeCourseData(data, range[0], range[1]);
      ui->widgetPlot->SetCurrentFrame(overlay->GetActiveFrame());
      setWindowTitle(QString("Time Course (%1)").arg(overlay->GetName()));
    }
  }
}

void WindowTimeCourse::OnFrameChanged(int frame)
{
  QString type = MainWindow::GetMainWindow()->GetCurrentLayerType();
  if (type != "MRI" && type != "Surface")
    type == "MRI";
  if (type == "MRI")
  {
    LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer("MRI"));
    if (layer && frame != layer->GetActiveFrame() && frame < layer->GetNumberOfFrames())
    {
      layer->SetActiveFrame(frame);
    }
  }
}

void WindowTimeCourse::SetCurrentFrame(int n)
{
  ui->widgetPlot->SetCurrentFrame(n);
}

void WindowTimeCourse::OnLayerCorrelationSurfaceChanged()
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer("MRI"));
  if (layer && layer->GetCorrelationSurface())
  {
    hide();
  }
}
