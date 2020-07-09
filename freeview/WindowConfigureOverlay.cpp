/*
 * Original Author: Ruopeng Wang
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
#include "WindowConfigureOverlay.h"
#include "ui_WindowConfigureOverlay.h"
#include "LayerSurface.h"
#include "SurfaceOverlayProperty.h"
#include "SurfaceOverlay.h"
#include "SurfaceLabel.h"
#include "LayerPropertySurface.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>
#include <QSettings>
#include <QTimer>
#include "DialogScreenshotOverlay.h"

WindowConfigureOverlay::WindowConfigureOverlay(QWidget *parent) :
  QWidget(parent), UIUpdateHelper(),
  ui(new Ui::WindowConfigureOverlay),
  m_dSavedOffset(0)
{
  ui->setupUi(this);
  ui->layoutOverlayList->removeWidget(ui->labelShortCut);
  setWindowFlags( Qt::Tool );
  m_fDataCache = NULL;
  ui->widgetHistogram->SetNumberOfBins( 200 );
  ui->widgetHolderAddPoint->hide();
  ui->checkBoxClearLower->hide();
  ui->checkBoxClearHigher->hide();
  ui->pushButtonFlip->hide();
  ui->pushButtonLoadCustom->hide();
  ui->pushButtonSaveCustom->hide();
  ui->widgetColorPicker->setCurrentColor(Qt::green);
  m_rangeOverall[0] = 0;
  m_rangeOverall[1] = 1;
  connect(ui->widgetHistogram, SIGNAL(MarkerChanged()), this, SLOT(OnHistogramMarkerChanged()));
  connect(ui->checkBoxAutoApply, SIGNAL(toggled(bool)), this, SLOT(CheckApply(bool)));
  connect(ui->checkBoxApplyToAll, SIGNAL(toggled(bool)), this, SLOT(CheckApply(bool)));
  connect(ui->checkBoxAutoFrame, SIGNAL(toggled(bool)), this, SLOT(OnCheckAutoFrameByVertex(bool)));
  connect(ui->pushButtonApply, SIGNAL(clicked(bool)), SLOT(OnButtonClicked()));
  connect(ui->pushButtonCancel, SIGNAL(clicked(bool)), SLOT(OnButtonClicked()));
  connect(ui->pushButtonScreenshot, SIGNAL(clicked(bool)), SLOT(OnButtonClicked()));
  connect(ui->pushButtonHelp, SIGNAL(clicked(bool)), SLOT(OnButtonClicked()));
  connect(ui->checkBoxFixedAxes, SIGNAL(toggled(bool)), SLOT(OnCheckFixedAxes(bool)));
  connect(ui->pushButtonLoadCustom, SIGNAL(clicked(bool)), SLOT(OnButtonLoadCustom()));
  connect(ui->pushButtonSaveCustom, SIGNAL(clicked(bool)), SLOT(OnButtonSaveCustom()));

  m_layerSurface = NULL;
  QSettings settings;
  QVariant v = settings.value("WindowConfigureOverlay/Geometry");
  if (v.isValid())
  {
    this->restoreGeometry(v.toByteArray());
  }
  v = settings.value("WindowConfigureOverlay/AutoApply");
  if (!v.isValid())
    v = true;
  ui->checkBoxAutoApply->setChecked(v.toBool());
  ui->checkBoxAutoFrame->setChecked(settings.value("WindowConfigureOverlay/AutoFrame").toBool());
  ui->checkBoxFixedAxes->setChecked(settings.value("WindowConfigureOverlay/FixedAxes", true).toBool());

  LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection("MRI");
  connect(lc, SIGNAL(LayerAdded(Layer*)), this, SLOT(UpdateUI()));
  connect(lc, SIGNAL(LayerRemoved(Layer*)), this, SLOT(UpdateUI()));
  connect(MainWindow::GetMainWindow(), SIGNAL(CycleOverlayRequested()), SLOT(OnCycleOverlay()));

  m_dlgScreenshot = new DialogScreenshotOverlay(this);
  m_dlgScreenshot->hide();
}

WindowConfigureOverlay::~WindowConfigureOverlay()
{
  if (m_fDataCache)
    delete[] m_fDataCache;
  m_fDataCache = 0;

  QSettings settings;
  settings.setValue("WindowConfigureOverlay/Geometry", this->saveGeometry());
  settings.setValue("WindowConfigureOverlay/AutoApply", ui->checkBoxAutoApply->isChecked());
  settings.setValue("WindowConfigureOverlay/AutoFrame", ui->checkBoxAutoFrame->isChecked());
  settings.setValue("WindowConfigureOverlay/FixedAxes", ui->checkBoxFixedAxes->isChecked());

  delete ui;
}

void WindowConfigureOverlay::showEvent(QShowEvent *)
{
  UpdateUI();
  UpdateGraph();
  UpdateGeometry();
}

void WindowConfigureOverlay::hideEvent(QHideEvent *)
{
  m_dlgScreenshot->hide();
}

void WindowConfigureOverlay::resizeEvent(QResizeEvent *e)
{
  UpdateGeometry();
}

void WindowConfigureOverlay::UpdateGeometry()
{
  QRect rc = ui->labelShortCut->geometry();
  rc.moveLeft(ui->comboBoxOverlayList->geometry().right()+10);
  rc.moveCenter(QPoint(rc.center().x(), ui->comboBoxOverlayList->geometry().center().y()+1));
  ui->labelShortCut->setGeometry(rc);
}

void WindowConfigureOverlay::OnActiveSurfaceChanged(Layer* layer)
{
  if (m_layerSurface)
  {
    disconnect(m_layerSurface, 0, this, 0);
  }
  m_layerSurface = qobject_cast<LayerSurface*>(layer);
  if (m_layerSurface)
  {
    disconnect(m_layerSurface, 0, this, 0);
    connect(m_layerSurface, SIGNAL(SurfaceOverlyDataUpdated()),
            this, SLOT(UpdateUI()), Qt::UniqueConnection);
    connect(m_layerSurface, SIGNAL(SurfaceOverlyDataUpdated()),
            this, SLOT(UpdateGraph()), Qt::QueuedConnection);
    connect(m_layerSurface, SIGNAL(ActiveOverlayChanged(int)),
            this, SLOT(OnActiveOverlayChanged()), Qt::UniqueConnection);
    connect(m_layerSurface, SIGNAL(SurfaceLabelAdded(SurfaceLabel*)),
            this, SLOT(OnSurfaceLabelAdded(SurfaceLabel*)), Qt::UniqueConnection);
    connect(m_layerSurface, SIGNAL(SurfaceLabelDeleted(SurfaceLabel*)),
            this, SLOT(UpdateUI()), Qt::UniqueConnection);
  }

  OnActiveOverlayChanged();
}

void WindowConfigureOverlay::OnActiveOverlayChanged()
{
  if (m_fDataCache)
    delete[] m_fDataCache;
  m_fDataCache = 0;

  UpdateUI();
  OnCheckFixedAxes(ui->checkBoxFixedAxes->isChecked(), false);
  UpdateGraph();
}

void WindowConfigureOverlay::UpdateUI()
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    QList<QWidget*> allwidgets = this->findChildren<QWidget*>();
    for ( int i = 0; i < allwidgets.size(); i++ )
    {
      allwidgets[i]->blockSignals( true );
    }

    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    SurfaceOverlayProperty* p = overlay->GetProperty();

    ui->comboBoxOverlayList->clear();
    for (int i = 0; i < m_layerSurface->GetNumberOfOverlays(); i++)
    {
      SurfaceOverlay* ol = m_layerSurface->GetOverlay(i);
      ui->comboBoxOverlayList->addItem(ol->GetName());
      if (ol == overlay)
        ui->comboBoxOverlayList->setCurrentIndex(i);
    }

    ui->sliderOpacity->setValue( (int)( p->GetOpacity() * 100 ) );
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, p->GetOpacity() );

    ui->sliderFrame->setRange(0, overlay->GetNumberOfFrames()-1);
    ui->sliderFrame->setValue(overlay->GetActiveFrame());
    ui->spinBoxFrame->setRange(0, overlay->GetNumberOfFrames()-1);
    ui->spinBoxFrame->setValue(overlay->GetActiveFrame());
    ui->labelFrameRange->setText(QString("0-%1").arg(overlay->GetNumberOfFrames()-1));
    ui->groupBoxFrame->setVisible(overlay->GetNumberOfFrames() > 1);

    //   ui->radioButtonGreenRed ->setChecked( p->GetColorScale() == SurfaceOverlayProperty::CS_GreenRed );
    //   ui->radioButtonBlueRed    ->setChecked( p->GetColorScale() == SurfaceOverlayProperty::CS_BlueRed );
    ui->radioButtonHeat       ->setChecked( p->GetColorScale() == SurfaceOverlayProperty::CS_Heat );
    ui->radioButtonColorWheel ->setChecked( p->GetColorScale() == SurfaceOverlayProperty::CS_ColorWheel );
    ui->radioButtonCustom  ->setChecked( p->GetColorScale() == SurfaceOverlayProperty::CS_Custom );

    ui->pushButtonLoadCustom->setVisible(ui->radioButtonCustom->isChecked());
    ui->pushButtonSaveCustom->setVisible(ui->radioButtonCustom->isChecked());
    ui->pushButtonFlip->setVisible(ui->radioButtonCustom->isChecked());
    ui->checkBoxClearHigher->setVisible(ui->radioButtonCustom->isChecked());
    ui->checkBoxClearLower->setVisible(ui->radioButtonCustom->isChecked());

    ui->checkBoxUsePercentile->setChecked(p->GetUsePercentile());
    ui->widgetHistogram->SetUsePercentile(p->GetUsePercentile());
    ui->checkBoxUseNonZeroVertices->setChecked(p->GetIgnoreZeros());
    ui->checkBoxUseNonZeroVertices->setVisible(p->GetUsePercentile());

    ui->checkBoxMaskInverse->setChecked(p->GetMaskInverse());

    if (p->GetUsePercentile())
    {
      ChangeLineEditNumber( ui->lineEditMin, overlay->PositionToPercentile(p->GetMinPoint()) );
      ChangeLineEditNumber( ui->lineEditMid, overlay->PositionToPercentile(p->GetMidPoint()) );
      ChangeLineEditNumber( ui->lineEditMax, overlay->PositionToPercentile(p->GetMaxPoint()) );
    }
    else
    {
      ChangeLineEditNumber( ui->lineEditMin, p->GetMinPoint() );
      ChangeLineEditNumber( ui->lineEditMid, p->GetMidPoint() );
      ChangeLineEditNumber( ui->lineEditMax, p->GetMaxPoint() );
    }
    ChangeLineEditNumber( ui->lineEditOffset, p->GetOffset());
    m_dSavedOffset = p->GetOffset();

    ui->radioButtonLinear        ->setChecked( p->GetColorMethod() == SurfaceOverlayProperty::CM_Linear );
    ui->radioButtonLinearOpaque  ->setChecked( p->GetColorMethod() == SurfaceOverlayProperty::CM_LinearOpaque );
    ui->radioButtonPiecewise     ->setChecked( p->GetColorMethod() == SurfaceOverlayProperty::CM_Piecewise );

    ui->checkBoxInverse->setChecked( p->GetColorInverse() );
    ui->checkBoxTruncate->setChecked( p->GetColorTruncate() );
    ui->checkBoxClearLower->setChecked(p->GetClearLower());
    ui->checkBoxClearHigher->setChecked(p->GetClearHigher());

    ui->labelMid->setEnabled( ui->radioButtonPiecewise->isChecked() );
    ui->lineEditMid->setEnabled( ui->radioButtonPiecewise->isChecked() );

    ui->checkBoxEnableSmooth->setChecked(p->GetSmooth());
    ui->spinBoxSmoothSteps->setValue(p->GetSmoothSteps());
    ui->spinBoxSmoothSteps->setEnabled(p->GetSmooth());

    QGradientStops stops = p->GetCustomColorScale();
    m_markers.clear();
    for (int i = 0; i < stops.size(); i++)
    {
      LineMarker m;
      m.position = stops[i].first;
      m.color = stops[i].second;
      m_markers << m;
    }

    int nFrames = overlay->GetNumberOfFrames();
    if (nFrames == 1)
      ui->groupBoxFrame->setEnabled(true);
    ui->checkBoxComputeCorrelation->setChecked(overlay->GetComputeCorrelation());

    ui->comboBoxVolumes->clear();
    ui->comboBoxVolumes->addItem("Self");
    if (!overlay->GetCorrelationSourceVolume())
      ui->comboBoxVolumes->setCurrentIndex(0);
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
    foreach (Layer* mri, layers)
    {
      if ((qobject_cast<LayerMRI*>(mri))->GetNumberOfFrames() == nFrames)
      {
        ui->comboBoxVolumes->addItem(mri->GetName(), QVariant::fromValue((QObject*)mri));
        if (overlay->GetCorrelationSourceVolume() == mri)
        {
          ui->comboBoxVolumes->setCurrentIndex(ui->comboBoxVolumes->count()-1);
        }
      }
    }

    ui->comboBoxMask->clear();
    ui->comboBoxMask->addItem("None");
    int nSel = 0;
    for (int i = 0; i < m_layerSurface->GetNumberOfLabels(); i++)
    {
      SurfaceLabel* label = m_layerSurface->GetLabel(i);
      ui->comboBoxMask->addItem(label->GetName(), QVariant::fromValue((QObject*)label));
      if (label == overlay->GetProperty()->GetMask())
        nSel = i+1;
    }
    ui->comboBoxMask->addItem("Load...");
    ui->comboBoxMask->setCurrentIndex(nSel);
    ui->checkBoxMaskInverse->setEnabled(nSel > 0);

    ui->checkBoxComputeCorrelation->setEnabled(ui->comboBoxVolumes->count() > 0);
    ui->groupBoxCorrelation->setVisible(nFrames > 1 && ui->comboBoxVolumes->count() > 0);

    for ( int i = 0; i < allwidgets.size(); i++ )
    {
      allwidgets[i]->blockSignals( false );
    }
  }
  else
  {
    close();
  }
}

void WindowConfigureOverlay::OnButtonClicked()
{
  if (sender() == ui->pushButtonHelp)
  {
    QMessageBox::information(this, "Help", "Drag the handle to move point.\n\nAt Custom mode:\nDouble-click on the handle to change point color.\nShift+Click on the handle to remove point.");
  }
  else if (sender() == ui->pushButtonApply)
  {
    OnApply();
  }
  else if (sender() == ui->pushButtonCancel)
  {
    close();
  }
  else if (sender() == ui->pushButtonScreenshot)
  {
    m_dlgScreenshot->show();
    m_dlgScreenshot->raise();
  }
}

void WindowConfigureOverlay::OnApply()
{
  if ( !m_layerSurface || !m_layerSurface->GetActiveOverlay() )
  {
    return;
  }

  SurfaceOverlayProperty* p = m_layerSurface->GetActiveOverlay()->GetProperty();
  bool smooth_changed = (p->GetSmooth() != ui->checkBoxEnableSmooth->isChecked() ||
      p->GetSmoothSteps() != ui->spinBoxSmoothSteps->value() );
  if ( UpdateOverlayProperty( p ) )
  {
    if (smooth_changed)
      m_layerSurface->GetActiveOverlay()->UpdateSmooth();
    else
      p->EmitColorMapChanged();
    if (ui->checkBoxApplyToAll->isChecked())
    {
      for (int i = 0; i < m_layerSurface->GetNumberOfOverlays(); i++)
      {
        SurfaceOverlay* so = m_layerSurface->GetOverlay(i);
        if (so != m_layerSurface->GetActiveOverlay())
        {
          smooth_changed = (so->GetProperty()->GetSmooth() != ui->checkBoxEnableSmooth->isChecked() ||
              so->GetProperty()->GetSmoothSteps() != ui->spinBoxSmoothSteps->value() );
          so->GetProperty()->Copy(p);
          if (smooth_changed)
            so->UpdateSmooth();
          else
            so->GetProperty()->EmitColorMapChanged();
        }
      }
    }
  }
}

void WindowConfigureOverlay::CheckApply(bool bChecked)
{
  if (bChecked)
  {
    OnApply();
  }
}

bool WindowConfigureOverlay::UpdateOverlayProperty( SurfaceOverlayProperty* p )
{
  if (!m_layerSurface || !m_layerSurface->GetActiveOverlay())
    return false;

  p->SetOpacity( ui->doubleSpinBoxOpacity->value() );
  p->SetUsePercentile(ui->checkBoxUsePercentile->isChecked());
  p->SetIgnoreZeros(ui->checkBoxUseNonZeroVertices->isChecked());

  bool bOK;
  double dValue = ui->lineEditMin->text().toDouble(&bOK);
  SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
  if ( bOK )
  {
    if (ui->checkBoxUsePercentile->isChecked())
      dValue = overlay->PercentileToPosition(dValue, p->GetIgnoreZeros());
    p->SetMinPoint( dValue );
  }
  else
  {
    return false;
  }

  dValue = ui->lineEditMid->text().toDouble(&bOK);
  if ( bOK )
  {
    if (ui->checkBoxUsePercentile->isChecked())
      dValue = overlay->PercentileToPosition(dValue, p->GetIgnoreZeros());
    p->SetMidPoint( dValue );
  }
  else
  {
    return false;
  }

  dValue = ui->lineEditMax->text().toDouble(&bOK);
  if ( bOK )
  {
    if (ui->checkBoxUsePercentile->isChecked())
      dValue = overlay->PercentileToPosition(dValue, p->GetIgnoreZeros());
    p->SetMaxPoint( dValue );
  }
  else
  {
    return false;
  }

  dValue = ui->lineEditOffset->text().toDouble(&bOK);
  if ( bOK )
  {
    p->SetOffset( dValue );
  }
  else
  {
    return false;
  }

  if ( ui->radioButtonHeat->isChecked() )
  {
    p->SetColorScale( SurfaceOverlayProperty::CS_Heat );
  }
  else if ( ui->radioButtonColorWheel->isChecked() )
  {
    p->SetColorScale( SurfaceOverlayProperty::CS_ColorWheel );
  }
  else if ( ui->radioButtonCustom->isChecked() )
  {
    p->SetColorScale( SurfaceOverlayProperty::CS_Custom );
  }

  if ( ui->radioButtonLinear->isChecked() )
  {
    p->SetColorMethod( SurfaceOverlayProperty::CM_Linear );
  }
  else if ( ui->radioButtonLinearOpaque->isChecked() )
  {
    p->SetColorMethod( SurfaceOverlayProperty::CM_LinearOpaque );
  }
  else if ( ui->radioButtonPiecewise->isChecked() )
  {
    p->SetColorMethod( SurfaceOverlayProperty::CM_Piecewise );
  }

  p->SetColorInverse( ui->checkBoxInverse->isChecked() );
  p->SetColorTruncate( ui->checkBoxTruncate->isChecked() );
  p->SetClearLower(ui->checkBoxClearLower->isChecked());
  p->SetClearHigher(ui->checkBoxClearHigher->isChecked());

  QGradientStops stops;
  for (int i = 0; i < m_markers.size(); i++)
  {
    stops << QGradientStop(m_markers[i].position, m_markers[i].color);
  }
  p->SetCustomColorScale(stops);

  p->SetSmooth(ui->checkBoxEnableSmooth->isChecked());
  p->SetSmoothSteps(ui->spinBoxSmoothSteps->value());

  p->SetMask(qobject_cast<SurfaceLabel*>(ui->comboBoxMask->itemData(ui->comboBoxMask->currentIndex()).value<QObject*>()));
  p->SetMaskInverse(ui->checkBoxMaskInverse->isChecked());

  return true;
}

void WindowConfigureOverlay::UpdateGraph(bool bApply)
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    if ( overlay )
    {
      double range[2];
      if (ui->checkBoxFixedAxes->isChecked())
      {
        range[0] = m_rangeOverall[0];
        range[1] = m_rangeOverall[1];
      }
      else
        overlay->GetDisplayRange( range );
      if (range[0] == range[1])
      {
        return;
      }

      SurfaceOverlayProperty* p = new SurfaceOverlayProperty( overlay );
      UpdateOverlayProperty( p );
      if (m_fDataCache)
        ui->widgetHistogram->SetInputData( m_fDataCache, overlay->GetDataSize(), range);
      else
        ui->widgetHistogram->SetInputData( overlay->GetData(), overlay->GetDataSize(), range);
      ui->widgetHistogram->SetSymmetricMarkers(p->GetColorScale() <= SurfaceOverlayProperty::CS_BlueRed);
      ui->widgetHistogram->SetMarkerEditable(p->GetColorScale() == SurfaceOverlayProperty::CS_Custom);

      int nBins = ui->widgetHistogram->GetNumberOfBins();
      float* fData = new float[ nBins ];
      unsigned char* nColorTable = new unsigned char[ nBins*4 ];

      ui->widgetHistogram->GetOutputRange(range);
      double bin_width = ( range[1] - range[0] ) / nBins;
      int rgb[3];
      double* dColor = m_layerSurface->GetProperty()->GetBinaryColor();
      rgb[0] = (int)( dColor[0] * 255 );
      rgb[1] = (int)( dColor[1] * 255 );
      rgb[2] = (int)( dColor[2] * 255 );
      for ( int i = 0; i < nBins; i++ )
      {
        nColorTable[i*4] = rgb[0];
        nColorTable[i*4+1] = rgb[1];
        nColorTable[i*4+2] = rgb[2];
        nColorTable[i*4+3] = 255;

        fData[i] = range[0] + ( i + 0.5 ) * bin_width;
      }
      p->MapOverlayColor( fData, nColorTable, nBins );
      ui->widgetHistogram->SetColorTableData( nColorTable, false );
      delete[] fData;
      delete[] nColorTable;

      LineMarkers markers;
      if (p->GetColorScale() == SurfaceOverlayProperty::CS_Custom)
      {
        markers = m_markers;
      }
      else
      {
        // rebuild marker lines for display
        LineMarker marker;
        marker.position = p->GetMinPoint()+p->GetOffset();
        marker.color = QColor( 255, 0, 0 );
        marker.movable = true;
        markers.push_back( marker );

        if ( p->GetColorMethod() == SurfaceOverlayProperty::CM_Piecewise )
        {
          marker.position = p->GetMidPoint()+p->GetOffset();
          marker.color = QColor( 0, 0, 255 );
          markers.push_back( marker );
        }

        marker.position = p->GetMaxPoint()+p->GetOffset();
        marker.color = QColor( 0, 215, 0 );
        markers.push_back( marker );
      }
      ui->widgetHistogram->SetMarkers( markers );
      delete p;
    }
    if (bApply && ui->checkBoxAutoApply->isChecked())
      OnApply();
  }
}

void WindowConfigureOverlay::OnSliderOpacity( int nVal )
{
  ui->doubleSpinBoxOpacity->blockSignals( true );
  ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, nVal / 100.0 );
  ui->doubleSpinBoxOpacity->blockSignals(false);
  UpdateGraph(true);
}

void WindowConfigureOverlay::OnSpinBoxOpacity( double dVal )
{
  ui->sliderOpacity->blockSignals(true);
  ui->sliderOpacity->setValue( (int)(dVal * 100) );
  ui->sliderOpacity->blockSignals(false);
  UpdateGraph(true);
}

void WindowConfigureOverlay::UpdateThresholdChanges()
{
  if ( !ui->radioButtonPiecewise->isChecked() )   // do not adjust mid point automatically in Piecewise mode
  {
    bool bOK;
    double dmin = ui->lineEditMin->text().trimmed().toDouble(&bOK);
    double dmax = ui->lineEditMax->text().trimmed().toDouble(&bOK);

    if ( bOK )
    {
      ui->lineEditMid->blockSignals(true);
      ChangeLineEditNumber(ui->lineEditMid, ( dmax + dmin ) / 2 );
      ui->lineEditMid->blockSignals(false);
    }
  }
  UpdateGraph(true);
}

void WindowConfigureOverlay::OnHistogramMouseButtonPressed(int button, double value)
{
  if (!m_layerSurface || !m_layerSurface->GetActiveOverlay())
    return;

  SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();

  value -= m_dSavedOffset;
  if (ui->checkBoxUsePercentile->isChecked())
    value = overlay->PositionToPercentile(value, ui->checkBoxUseNonZeroVertices->isChecked());

  if (button == Qt::LeftButton)
  {
    ChangeLineEditNumber(ui->lineEditMin, value);
  }
  else if (button == Qt::MidButton)
  {
    ChangeLineEditNumber(ui->lineEditMid, value);
  }
  else if (button == Qt::RightButton)
  {
    ChangeLineEditNumber(ui->lineEditMax, value);
  }
  UpdateThresholdChanges();
}

void WindowConfigureOverlay::OnHistogramMarkerChanged()
{
  if (!m_layerSurface || !m_layerSurface->GetActiveOverlay())
    return;

  SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();

  LineMarkers markers = ui->widgetHistogram->GetMarkers();
  if (ui->radioButtonCustom->isChecked())
  {
    m_markers = markers;
    UpdateGraph(true);
  }
  else
  {
    bool bUsePercentile = ui->checkBoxUsePercentile->isChecked();
    bool bIgnoreZeros = ui->checkBoxUseNonZeroVertices->isChecked();
    for (int i = 0; i < markers.size(); i++)
    {
      if (i == 0)
      {
        if (bUsePercentile)
          ChangeLineEditNumber(ui->lineEditMin, overlay->PositionToPercentile(markers[i].position-m_dSavedOffset, bIgnoreZeros),
                               2, true);
        else
          ChangeLineEditNumber(ui->lineEditMin, markers[i].position-m_dSavedOffset, 2, true);
      }
      if (i == 1)
      {
        if (ui->radioButtonPiecewise->isChecked() && ui->radioButtonHeat->isChecked())
        {
          if (bUsePercentile)
            ChangeLineEditNumber(ui->lineEditMid, overlay->PositionToPercentile(markers[i].position-m_dSavedOffset, bIgnoreZeros),
                                 2, true);
          else
            ChangeLineEditNumber(ui->lineEditMid, markers[i].position-m_dSavedOffset, 2, true);
        }
      }
    }
    if (bUsePercentile)
      ChangeLineEditNumber(ui->lineEditMax, overlay->PositionToPercentile(markers[markers.size()-1].position-m_dSavedOffset, bIgnoreZeros),
          2, true);
    else
      ChangeLineEditNumber(ui->lineEditMax, markers[markers.size()-1].position-m_dSavedOffset, 2, true);
    UpdateThresholdChanges();
  }
}

void WindowConfigureOverlay::OnButtonAdd()
{
  bool bOK;
  double pos = ui->lineEditNewPoint->text().toDouble(&bOK);
  if (!bOK)
  {
    QMessageBox::warning(this, "Error", "Please enter valid new point position.");
    return;
  }
  double range[2];
  ui->widgetHistogram->GetOutputRange(range);
  if (pos < range[0] || pos > range[1])
  {
    if (pos < range[0])
      range[0] = pos;
    else
      range[1] = pos;
    if (m_layerSurface)
    {
      SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
      if (overlay)
        overlay->SetDisplayRange(range);
      OnCheckFixedAxes(ui->checkBoxFixedAxes->isChecked(), false);
    }
    //    QMessageBox::warning(this, "Error", "New point out of range.");
    //    return;
  }
  ui->widgetHistogram->AddMarker(pos, ui->widgetColorPicker->currentColor());
}

void WindowConfigureOverlay::OnSmoothChanged()
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    if ( overlay )
    {
      if (!m_fDataCache)
        m_fDataCache = new float[overlay->GetDataSize()];

      if (ui->checkBoxEnableSmooth->isChecked())
      {
        overlay->SmoothData(ui->spinBoxSmoothSteps->value(), m_fDataCache);
      }
      else
      {
        memcpy(m_fDataCache, overlay->GetUnsmoothedData(), sizeof(float)*overlay->GetDataSize());
      }
      UpdateGraph(true);
    }
  }
}

void WindowConfigureOverlay::OnTextThresholdChanged(const QString &strg)
{
  if (!m_layerSurface || !m_layerSurface->GetActiveOverlay())
    return;

  SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();

  bool ok;
  double val = strg.toDouble(&ok);
  if (!ok)
    return;

  bool bIgnoreZeros = ui->checkBoxUseNonZeroVertices->isChecked();
  double dOffset = ui->lineEditOffset->text().trimmed().toDouble(&ok);
  if (!ok)
    dOffset = m_dSavedOffset;

  LineMarkers markers = ui->widgetHistogram->GetMarkers();
  if (markers.isEmpty())
    return;

  if (sender() == ui->lineEditMax || sender() == NULL)
  {
    LineMarker marker = markers.last();
    if (ui->checkBoxUsePercentile->isChecked())
      marker.position = overlay->PercentileToPosition(val, bIgnoreZeros)+dOffset;
    else
      marker.position = val+dOffset;
    markers[markers.size()-1] = marker;
    double val2 = ui->lineEditMin->text().toDouble(&ok);
    if (markers.size() == 2 && ok)
    {
      this->ChangeLineEditNumber(ui->lineEditMid, (val+val2)/2);
    }
  }
  else if (sender() == ui->lineEditMin)
  {
    LineMarker marker = markers.first();
    if (ui->checkBoxUsePercentile->isChecked())
      marker.position = overlay->PercentileToPosition(val, bIgnoreZeros)+dOffset;
    else
      marker.position = val+dOffset;
    markers[0] = marker;
    double val2 = ui->lineEditMax->text().toDouble(&ok);
    if (markers.size() == 2 && ok)
    {
      this->ChangeLineEditNumber(ui->lineEditMid, (val+val2)/2);
    }
  }
  else if (sender() == ui->lineEditOffset)
  {
    for (int i = 0; i < markers.size(); i++)
    {
      LineMarker marker = markers[i];
      marker.position += (val - m_dSavedOffset);
      markers[i] = marker;
    }
    m_dSavedOffset = val;
  }
  ui->widgetHistogram->SetMarkers(markers);
  UpdateGraph(true);
}

void WindowConfigureOverlay::OnFrameChanged(int nFrame)
{
  if (sender() != ui->spinBoxFrame)
  {
    ui->spinBoxFrame->blockSignals(true);
    ui->spinBoxFrame->setValue(nFrame);
    ui->spinBoxFrame->blockSignals(false);
  }
  if (sender() != ui->sliderFrame)
  {
    ui->sliderFrame->blockSignals(true);
    ui->sliderFrame->setValue(nFrame);
    ui->sliderFrame->blockSignals(false);
  }
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    overlay->SetActiveFrame(nFrame);
    delete[] m_fDataCache;
    m_fDataCache = NULL;
    UpdateGraph(true);
    emit ActiveFrameChanged(nFrame);
  }
}

void WindowConfigureOverlay::OnCurrentVertexChanged()
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay()
       && m_layerSurface->GetActiveOverlay()->GetNumberOfFrames() > 1 )
  {
    int nVertex;
    if (m_layerSurface->IsInflated())
      nVertex = m_layerSurface->GetCurrentVertex();
    else
      nVertex = m_layerSurface->GetVertexIndexAtTarget( m_layerSurface->GetSlicePosition(), NULL );

    if (nVertex >= 0 && ui->checkBoxAutoFrame->isChecked()
        && nVertex < m_layerSurface->GetActiveOverlay()->GetNumberOfFrames())
    {
      OnFrameChanged(nVertex);

      // even if not auto apply, still call apply
      if (!ui->checkBoxAutoApply->isChecked())
        OnApply();
    }
  }
}

void WindowConfigureOverlay::OnCheckComputeCorrelation(bool bChecked)
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    if (bChecked)
    {
      int n = ui->comboBoxVolumes->currentIndex();
      if (n >= 0)
      {
        OnComboCorrelationVolume(n);
      }
    }
    overlay->SetComputeCorrelation(bChecked);
    delete[] m_fDataCache;
    m_fDataCache = NULL;

    UpdateGraph();
  }
}

void WindowConfigureOverlay::OnComboCorrelationVolume(int n)
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    LayerMRI* mri = qobject_cast<LayerMRI*>(ui->comboBoxVolumes->itemData(n).value<QObject*>());
    overlay->SetCorrelationSourceVolume(mri);
  }
}

void WindowConfigureOverlay::OnCheckUsePercentile(bool bChecked)
{
  ui->checkBoxUseNonZeroVertices->setVisible(bChecked);
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    overlay->GetProperty()->SetUsePercentile(bChecked);
    ui->widgetHistogram->SetUsePercentile(bChecked);
    UpdateUI();
  }
}

void WindowConfigureOverlay::OnCheckUseNonZeroVertices(bool bChecked)
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    overlay->GetProperty()->SetIgnoreZeros(bChecked);
    OnHistogramMarkerChanged();
    OnTextThresholdChanged(ui->lineEditMax->text());
  }
}

void WindowConfigureOverlay::OnCustomColorScale()
{
  ui->lineEditOffset->setText("0");
}

void WindowConfigureOverlay::OnComboMask(int n)
{
  if (n == ui->comboBoxMask->count()-1)
  {
    QString filename = QFileDialog::getOpenFileName( this, "Select label file",
                                                     MainWindow::GetMainWindow()->AutoSelectLastDir( "label" ),
                                                     "Label files (*)");
    if ( !filename.isEmpty())
    {
      LoadLabelMask(filename);
    }
    else
    {
      UpdateUI();
    }
  }
  else
  {
    if (ui->checkBoxAutoApply->isChecked())
      OnApply();
    UpdateUI();
  }
}

void WindowConfigureOverlay::LoadLabelMask(const QString& fn)
{
  setProperty("wait_for_label", true);
  emit MaskLoadRequested(fn);
}

void WindowConfigureOverlay::OnCheckInverseMask(bool bChecked)
{
  Q_UNUSED(bChecked);

  if (ui->checkBoxAutoApply->isChecked())
    OnApply();
}

void WindowConfigureOverlay::OnSurfaceLabelAdded(SurfaceLabel* label)
{
  Q_UNUSED(label);

  UpdateUI();
  if (property("wait_for_label").toBool())
  {
    setProperty("wait_for_label", false);
    ui->comboBoxMask->setCurrentIndex(1);
  }
}

void WindowConfigureOverlay::OnCheckAutoFrameByVertex(bool bChecked)
{
  if (bChecked)
    OnCurrentVertexChanged();
}

void WindowConfigureOverlay::OnComboOverlayChanged(int n)
{
  if (m_layerSurface && n < m_layerSurface->GetNumberOfOverlays())
  {
    m_layerSurface->SetActiveOverlay(n);
    emit OverlayChanged();
  }
}

void WindowConfigureOverlay::OnCycleOverlay()
{
  if (m_layerSurface && m_layerSurface->GetNumberOfOverlays() > 1)
  {
    ui->comboBoxOverlayList->setCurrentIndex((m_layerSurface->GetActiveOverlayIndex()+1)%m_layerSurface->GetNumberOfOverlays());
  }
}

void WindowConfigureOverlay::OnCheckFixedAxes(bool bChecked, bool bUpdateGraph)
{
  if (bChecked && m_layerSurface)
  {
    m_rangeOverall[0] = 1e10;
    m_rangeOverall[1] = -1e10;
    for (int i = 0; i < m_layerSurface->GetNumberOfOverlays(); i++)
    {
      SurfaceOverlay* ol = m_layerSurface->GetOverlay(i);
      double range[2];
      ol->GetDisplayRange(range);
      if (range[0] < m_rangeOverall[0])
        m_rangeOverall[0] = range[0];
      if (range[1] > m_rangeOverall[1])
        m_rangeOverall[1] = range[1];
    }
  }
  if (bUpdateGraph)
    UpdateGraph();
}

void WindowConfigureOverlay::OnButtonSaveCustom()
{
  QString filename = QFileDialog::getSaveFileName( this, "Save Color Scale",
                                                   MainWindow::GetMainWindow()->AutoSelectLastDir( "surf" ),
                                                   "All files (*)");
  if ( !filename.isEmpty() && m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    if (!overlay->GetProperty()->SaveCustomColorScale(filename))
      QMessageBox::warning(this, "Error", "Failed to save color scale to " + filename);
  }
}

void WindowConfigureOverlay::OnButtonLoadCustom()
{
  QString filename = QFileDialog::getOpenFileName( this, "Load Color Scale",
                                                   MainWindow::GetMainWindow()->AutoSelectLastDir( "surf" ),
                                                   "All files (*)");
  if ( !filename.isEmpty() && m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    if (!overlay->GetProperty()->LoadCustomColorScale(filename))
      QMessageBox::warning(this, "Error", "Failed to load color scale from " + filename);
    else
    {
      m_layerSurface->UpdateOverlay(true);
      overlay->EmitDataUpdated();
    }
  }
}
