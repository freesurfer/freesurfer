/**
 * @brief Tool window to display time course data
 *
 */
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
#include "FlowLayout.h"
#include <QColorDialog>
#include <QPointer>

ClickableLabel::ClickableLabel(QWidget* parent, Qt::WindowFlags f)
    : QLabel(parent)
{
}

ClickableLabel::~ClickableLabel() {}

void ClickableLabel::mousePressEvent(QMouseEvent* event)
{
    emit clicked();
}

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
  ui->widgetPlot->SetDarkMode(true);
  layoutLegend = new FlowLayout;
  ui->widgetLegend->setLayout(layoutLegend);

  connect(ui->widgetPlot, SIGNAL(FrameChanged(int)), this, SLOT(OnFrameChanged(int)), Qt::QueuedConnection);
  connect(ui->widgetPlot, SIGNAL(PlotRangeChanged()), this, SLOT(UpdateScaleInfo()), Qt::QueuedConnection);

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
}

void WindowTimeCourse::Clear()
{
  ui->widgetPlot->Clear();
  QLayoutItem* item;
  while ( ( item = layoutLegend->takeAt( 0 ) ) != NULL )
  {
    delete item->widget();
    delete item;
  }
}

void WindowTimeCourse::UpdateData(bool bForce)
{
  if (!isVisible() && !bForce)
    return;

  QString type = MainWindow::GetMainWindow()->GetCurrentLayerType();
  QList<QColor> colors;
  colors << Qt::yellow << Qt::cyan << Qt::red << Qt::magenta;
  Clear();
  if (type == "MRI")
  {
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
    for (int nl = layers.size()-1; nl >= 0; nl--)
    {
      LayerMRI* layer = qobject_cast<LayerMRI*>(layers[nl]);
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
        double val_min = vol->GetMinValue(), val_max = vol->GetFullMaxValue();
        QVariantMap info = layer->GetTimeSeriesInfo();
        TimeCourseData td;
        td.m_points = data;
        td.m_dMin = val_min;
        td.m_dMax = val_max;
        if (!layer->property("legend_color").value<QColor>().isValid())
          layer->setProperty("legend_color", colors[(layers.size()-1-nl)%colors.size()]);
        td.m_color = layer->property("legend_color").value<QColor>();
        td.m_strXUnit = info["unit"].toString();
        td.m_dXInterval = info["tr"].toDouble();
        td.m_dXOffset = info["offset"].toDouble();
        td.m_nId = layer->GetID();
        td.m_strName = layer->GetName();
        if (layer->property("timecourse_visible").isValid())
          td.m_bShow = layer->property("timecourse_visible").toBool();
        ui->widgetPlot->AddTimeCourseData(td);
        ui->widgetPlot->SetCurrentFrame(layer->GetActiveFrame());

        connect(layer, SIGNAL(CorrelationSurfaceChanged(LayerSurface*)),
                this, SLOT(OnLayerCorrelationSurfaceChanged()), Qt::UniqueConnection);

        QWidget* w = MakeLegendWidget(layer, td);
        layoutLegend->addWidget(w);
      }
      //  setWindowTitle(QString("Time Course (%1)").arg(layer->GetName()));
    }
  }
  else if (type == "Surface")
  {
    LayerSurface* surf = qobject_cast<LayerSurface*>(MainWindow::GetMainWindow()->GetActiveLayer(type));
    if (surf)
    {
      Clear();
      double pos[3];
      MainWindow::GetMainWindow()->GetLayerCollection("Surface")->GetSlicePosition(pos);
      int nVert = -1;
      if (surf->GetFileName().contains("inflated", Qt::CaseInsensitive) && surf->GetSourceSurface()->IsSurfaceLoaded(FSSurface::SurfaceWhite))
        nVert = surf->GetCurrentVertex();
      else
        nVert = surf->GetVertexIndexAtTarget(pos, NULL);
      if (nVert < 0)
        return;

      QList<SurfaceOverlay*> overlays = surf->GetOverlays();
      for (int no = 0; no < overlays.size(); no++)
      {
        if (overlays[no]->GetNumberOfFrames() <= 1)
        {
          overlays.removeAt(no);
          no--;
        }
      }

      for (int no = 0; no < overlays.size(); no++)
      {
        SurfaceOverlay* overlay = overlays[no];
        int nFrames = overlay->GetNumberOfFrames();
        float* buffer = new float[nFrames];
        double range[2];
        overlay->GetDataAtVertex(nVert, buffer);
        overlay->GetRawRange(range);
        QList<double> data;
        for (int i = 0; i < nFrames; i++)
          data << buffer[i];
        delete[] buffer;
        TimeCourseData td;
        td.m_points = data;
        td.m_dMin = range[0];
        td.m_dMax = range[1];
        if (!overlay->property("legend_color").value<QColor>().isValid())
          overlay->setProperty("legend_color", colors[(overlays.size()-1-no)%colors.size()]);
        if (overlay->property("timecourse_visible").isValid())
          td.m_bShow = overlay->property("timecourse_visible").toBool();
        td.m_color = overlay->property("legend_color").value<QColor>();
        td.m_nId = overlay->GetID();
        td.m_strName = overlay->GetName();
        ui->widgetPlot->AddTimeCourseData(td);
        ui->widgetPlot->SetCurrentFrame(overlay->GetActiveFrame());

        QWidget* w = MakeLegendWidget(overlay, td);
        layoutLegend->addWidget(w);
      }
      //  setWindowTitle(QString("Time Course (%1)").arg(overlay->GetName()));
    }
  }
}

QWidget* WindowTimeCourse::MakeLegendWidget(QObject* obj, const TimeCourseData& td)
{
  QWidget* w = new QWidget(this);
  QHBoxLayout* hbox = new QHBoxLayout;
  w->setLayout(hbox);
  QCheckBox* checkbox = new QCheckBox();
  checkbox->setChecked(td.m_bShow);
  checkbox->setProperty("data_id", td.m_nId);
  checkbox->setProperty("data_obj", qVariantFromValue(obj));
  connect(checkbox, SIGNAL(toggled(bool)), SLOT(OnCheckBoxShowData(bool)));
  checkbox->setCursor(Qt::PointingHandCursor);
  hbox->addWidget(checkbox);
  ClickableLabel* label = new ClickableLabel();
  label->setText(td.m_strName);
  label->setProperty("data_id", td.m_nId);
  hbox->addWidget(label);
  label->setStyleSheet(QString("color:rgb(%1,%2,%3)")
                          .arg(td.m_color.red()).arg(td.m_color.green()).arg(td.m_color.blue()));
  label->setProperty("data_obj", qVariantFromValue(obj));
  label->setCursor(Qt::PointingHandCursor);
  connect(label, SIGNAL(clicked()), SLOT(OnLegendLabelClicked()));

  return w;
}

void WindowTimeCourse::OnFrameChanged(int frame)
{
  QString type = MainWindow::GetMainWindow()->GetCurrentLayerType();
  if (type != "MRI" && type != "Surface")
    type = "MRI";
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

void WindowTimeCourse::OnCheckBoxShowData(bool bShow)
{
  QCheckBox* cb = qobject_cast<QCheckBox*>(sender());
  if (cb)
  {
    ui->widgetPlot->SetDataVisible(cb->property("data_id").toLongLong(), bShow);
    QObject* obj = cb->property("data_obj").value<QObject*>();
    if (obj)
      obj->setProperty("timecourse_visible", bShow);
  }
}

void WindowTimeCourse::OnLegendLabelClicked()
{
  QPointer<QLabel> l = qobject_cast<QLabel*>(sender());
  if (l)
  {
    QObject* obj = l->property("data_obj").value<QObject*>();
    QColor c = QColorDialog::getColor(obj?obj->property("legend_color").value<QColor>():Qt::white, this);
    if (c.isValid() && l)
    {
      l->setStyleSheet(QString("color:rgb(%1,%2,%3)")
                              .arg(c.red()).arg(c.green()).arg(c.blue()));
      obj->setProperty("legend_color", c);
      ui->widgetPlot->SetDataColor(l->property("data_id").toLongLong(), c);
    }
  }
}
