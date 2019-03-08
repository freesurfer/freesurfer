/**
 * @file  WindowTimeCourse.cpp
 * @brief Tool window to display time course data
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2014/04/29 18:10:42 $
 *    $Revision: 1.6 $
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
#include "SurfaceOverlayProperty.h"
#include "FSSurface.h"
#include <QSettings>
#include <QDebug>

WindowTimeCourse::WindowTimeCourse(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::WindowTimeCourse),
  lastMRI(NULL),
  lastSurface(NULL),
  lastOverlay(NULL)
{
  ui->setupUi(this);
  this->setWindowFlags(Qt::Tool);
  this->setWindowTitle("Time Course");
  ui->frameSecondPlot->hide();
  connect(ui->widgetPlot, SIGNAL(FrameChanged(int)), this, SLOT(OnFrameChanged(int)));
  connect(ui->widgetPlot, SIGNAL(PlotRangeChanged()), this, SLOT(UpdateScaleInfo()), Qt::QueuedConnection);
  connect(ui->comboBoxSecondPlot, SIGNAL(currentIndexChanged(int)), SLOT(OnComboSecondPlot(int)));

  QSettings s;
  QVariant v = s.value("WindowTimeCourse/Geomerty");
  if (v.isValid())
    this->restoreGeometry(v.toByteArray());
  ui->checkBoxAutoScale->setChecked(s.value("WindowTimeCourse/AutoScale", true).toBool());
  if (!ui->checkBoxAutoScale->isChecked())
    ui->checkBoxMaxScale->setChecked(true);
  ui->lineEditScale->setEnabled(false);
}

WindowTimeCourse::~WindowTimeCourse()
{
  QSettings s;
  s.setValue("WindowTimeCourse/Geomerty", this->saveGeometry());
  s.setValue("WindowTimeCourse/AutoScale", ui->checkBoxAutoScale->isChecked());
  delete ui;
}

void WindowTimeCourse::showEvent(QShowEvent *e)
{
  Q_UNUSED(e);
  UpdateData(true);
}

void WindowTimeCourse::UpdateUI()
{
  QString type = MainWindow::GetMainWindow()->GetCurrentLayerType();
  if (type == "MRI")
  {
    QList<LayerMRI*> valid_layers;
    LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer(type));
    if (layer && layer->GetNumberOfFrames() == 1 && lastMRI && MainWindow::GetMainWindow()->GetLayerCollection("MRI")->Contains(lastMRI))
      layer = lastMRI;
    if (layer)
    {
      QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
      LayerMRI* current_2nd_layer = qobject_cast<LayerMRI*>(
            ui->comboBoxSecondPlot->itemData(ui->comboBoxSecondPlot->currentIndex()).value<QObject*>());
      foreach (Layer* l, layers)
      {
        if (((LayerMRI*)l)->GetNumberOfFrames() == layer->GetNumberOfFrames() && l != layer)
          valid_layers << ((LayerMRI*)l);
      }
      ui->comboBoxSecondPlot->blockSignals(true);
      ui->comboBoxSecondPlot->clear();
      ui->comboBoxSecondPlot->addItem("Off");
      int n = 0;
      for (int i = 0; i < valid_layers.size(); i++)
      {
        if (valid_layers[i] == current_2nd_layer)
          n = i+1;
        ui->comboBoxSecondPlot->addItem(valid_layers[i]->GetName(), QVariant::fromValue((QObject*)valid_layers[i]));
      }
      ui->comboBoxSecondPlot->setCurrentIndex(n);
      ui->comboBoxSecondPlot->blockSignals(false);
    }
    ui->frameSecondPlot->setVisible(!valid_layers.isEmpty());
    ui->labelSecondPlot->setText("Second Volume");
  }
  else if (type == "Surface")
  {
    QList<SurfaceOverlay*> valid_overlays;
    LayerSurface* layer = qobject_cast<LayerSurface*>(MainWindow::GetMainWindow()->GetActiveLayer(type));
    if (layer && layer->GetNumberOfOverlays() > 1 && layer->GetActiveOverlay())
    {
      QList<SurfaceOverlay*> overlays = layer->GetOverlays();
      SurfaceOverlay* current_overlay = layer->GetActiveOverlay();
      SurfaceOverlay* current_2nd_overlay = qobject_cast<SurfaceOverlay*>(
            ui->comboBoxSecondPlot->itemData(ui->comboBoxSecondPlot->currentIndex()).value<QObject*>());
      foreach (SurfaceOverlay* l, overlays)
      {
        if (l->GetNumberOfFrames() == current_overlay->GetNumberOfFrames() && l != current_overlay)
          valid_overlays << l;
      }
      ui->comboBoxSecondPlot->blockSignals(true);
      ui->comboBoxSecondPlot->clear();
      ui->comboBoxSecondPlot->addItem("Off");
      int n = 0;
      for (int i = 0; i < valid_overlays.size(); i++)
      {
        if (valid_overlays[i] == current_2nd_overlay)
          n = i+1;
        ui->comboBoxSecondPlot->addItem(valid_overlays[i]->GetName(), QVariant::fromValue((QObject*)valid_overlays[i]));
      }
      ui->comboBoxSecondPlot->setCurrentIndex(n);
      ui->comboBoxSecondPlot->blockSignals(false);
    }
    ui->frameSecondPlot->setVisible(!valid_overlays.isEmpty());
    ui->labelSecondPlot->setText("Second Overlay");
  }
  else
    ui->comboBoxSecondPlot->hide();
}

void WindowTimeCourse::UpdateData(bool bForce)
{
  if (!isVisible() && !bForce)
    return;

  QString type = MainWindow::GetMainWindow()->GetCurrentLayerType();
  if (type == "MRI")
  {
    LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer(type));
    if (lastMRI == NULL)
    {
      QList<Layer*> mri_col = MainWindow::GetMainWindow()->GetLayers("MRI");
      foreach (Layer* mri, mri_col)
      {
        if (((LayerMRI*)mri)->GetNumberOfFrames() > 1)
        {
          lastMRI = ((LayerMRI*)mri);
          break;
        }
      }
    }
    if (layer && layer->GetNumberOfFrames() == 1 &&
        MainWindow::GetMainWindow()->GetLayerCollection("MRI")->Contains(lastMRI))
      layer = lastMRI;
    if (layer && layer->GetNumberOfFrames() > 1)
    {
      double ras[3];
      int n[3];
      MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetSlicePosition(ras);
      layer->RemapPositionToRealRAS(ras, ras);
      layer->RASToOriginalIndex(ras, n);
      QList<double> data, data2;
      for (int i = 0; i < layer->GetNumberOfFrames(); i++)
        data <<  layer->GetVoxelValueByOriginalIndex(n[0], n[1], n[2], i);
      LayerMRI* layer2 = qobject_cast<LayerMRI*>(ui->comboBoxSecondPlot->itemData(ui->comboBoxSecondPlot->currentIndex()).value<QObject*>());
      if (layer2)
      {
        MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetSlicePosition(ras);
        layer2->RemapPositionToRealRAS(ras, ras);
        layer2->RASToOriginalIndex(ras, n);
        for (int i = 0; i < layer2->GetNumberOfFrames(); i++)
        {
          data2 <<  layer2->GetVoxelValueByOriginalIndex(n[0], n[1], n[2], i);
        }
      }
      FSVolume* vol = layer->GetSourceVolume();
      double val_min = vol->GetMinValue(), val_max = vol->GetFullMaxValue();
      if (layer2)
      {
        val_min = qMin(val_min, layer2->GetSourceVolume()->GetMinValue());
        val_max = qMax(val_max, layer2->GetSourceVolume()->GetFullMaxValue());
      }
      ui->widgetPlot->SetTimeCourseData(data, val_min, val_max, layer->GetTR());
      ui->widgetPlot->SetSecondData(data2);
      ui->widgetPlot->SetCurrentFrame(layer->GetActiveFrame());

      QVariantMap info = layer->GetTimeSeriesInfo();
      ui->widgetPlot->SetXUnitInfo(info["tr"].toDouble(), info["offset"].toDouble(), info["unit"].toString());

      connect(layer, SIGNAL(CorrelationSurfaceChanged(LayerSurface*)),
              this, SLOT(OnLayerCorrelationSurfaceChanged()), Qt::UniqueConnection);
      setWindowTitle(QString("Time Course (%1)").arg(layer->GetName()));
      lastMRI = layer;
    }
  }
  else if (type == "Surface")
  {
    LayerSurface* surf = qobject_cast<LayerSurface*>(MainWindow::GetMainWindow()->GetActiveLayer(type));
    SurfaceOverlay* overlay = (surf ? surf->GetActiveOverlay() : NULL);
    if (surf && lastOverlay == NULL)
    {
      for (int i = 0; i < surf->GetNumberOfOverlays(); i++)
      {
        SurfaceOverlay* so = surf->GetOverlay(i);
        if (so->GetNumberOfFrames() > 1)
        {
          lastOverlay = so;
          lastSurface = surf;
          break;
        }
      }
    }
    if (overlay && overlay->GetNumberOfFrames() == 1 &&
        MainWindow::GetMainWindow()->GetLayerCollection("Surface")->Contains(lastSurface))
    {
      if (lastOverlay && lastOverlay->GetNumberOfFrames() > 1)
      {
        overlay = lastOverlay;
        surf = lastSurface;
      }
    }

    if (surf && overlay && overlay->GetNumberOfFrames() > 1)
    {
      double pos[3];
      MainWindow::GetMainWindow()->GetLayerCollection("Surface")->GetSlicePosition(pos);
      int nVert = -1;
      if (surf->GetFileName().contains("inflated", Qt::CaseInsensitive) && surf->GetSourceSurface()->IsSurfaceLoaded(FSSurface::SurfaceWhite))
        nVert = surf->GetCurrentVertex();
      else
        nVert = surf->GetVertexIndexAtTarget(pos, NULL);
      if (nVert < 0)
        return;

      int nFrames = overlay->GetNumberOfFrames();
      float* buffer = new float[nFrames];
      double range[2];
      overlay->GetDataAtVertex(nVert, buffer);
      overlay->GetRawRange(range);
      QList<double> data, data2;
      for (int i = 0; i < nFrames; i++)
        data << buffer[i];
      SurfaceOverlay* overlay2 = qobject_cast<SurfaceOverlay*>(ui->comboBoxSecondPlot->itemData(ui->comboBoxSecondPlot->currentIndex()).value<QObject*>());

      if (overlay2)
      {
        double range2[2];
        overlay2->GetDataAtVertex(nVert, buffer);
        overlay2->GetRawRange(range2);
        range[0] = qMin(range[0], range2[0]);
        range[1] = qMax(range[1], range2[1]);
        for (int i = 0; i < nFrames; i++)
          data2 << buffer[i];
      }
      delete[] buffer;
      ui->widgetPlot->SetTimeCourseData(data, range[0], range[1]);
      ui->widgetPlot->SetSecondData(data2);
      ui->widgetPlot->SetCurrentFrame(overlay->GetActiveFrame());
      setWindowTitle(QString("Time Course (%1)").arg(overlay->GetName()));
      lastSurface = surf;
      lastOverlay = overlay;
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
  else
  {
    emit OverlayFrameChanged(frame);
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

void WindowTimeCourse::UpdateScaleInfo()
{
  double range[2];
  ui->widgetPlot->GetPlotRange(range);
  ui->lineEditScale->blockSignals(true);
  ui->lineEditScale->setText(QString("%1, %2").arg(range[0]).arg(range[1]));
  ui->lineEditScale->blockSignals(false);
}

void WindowTimeCourse::OnCheckAutoScale(bool bChecked)
{
  ui->widgetPlot->SetAutoScale(bChecked);
  if (bChecked)
  {
    ui->checkBoxMaxScale->setChecked(false);
    ui->lineEditScale->setEnabled(false);
  }
  else if (!ui->checkBoxMaxScale->isChecked())
    ui->lineEditScale->setEnabled(true);
}

void WindowTimeCourse::OnCheckMaxScale(bool bChecked)
{
  if (bChecked)
  {
    ui->widgetPlot->ResetPlotRange();
    ui->checkBoxAutoScale->setChecked(false);
    ui->lineEditScale->setEnabled(false);
  }
  else if (!ui->checkBoxAutoScale->isChecked())
    ui->lineEditScale->setEnabled(true);
}

void WindowTimeCourse::OnLineEditScaleReturnPressed()
{
  QStringList list = ui->lineEditScale->text().trimmed().split(",", QString::SkipEmptyParts);
  if (list.size() != 2)
    list = ui->lineEditScale->text().trimmed().split(" ", QString::SkipEmptyParts);
  if (list.size() != 2)
    return;

  bool bOK;
  double range[2];
  range[0] = list[0].trimmed().toDouble(&bOK);
  if (!bOK)
    return;
  range[1] = list[1].trimmed().toDouble(&bOK);
  if (!bOK)
    return;

  if (range[1] <= range[0])
    range[1] = range[0]+1;
  ui->widgetPlot->SetPlotRange(range);
}

void WindowTimeCourse::OnComboSecondPlot(int nSel)
{
  UpdateData();
}

void WindowTimeCourse::OnCheckShowFrameNumber(bool b)
{
  ui->widgetPlot->SetShowFrameNumber(b);
}

