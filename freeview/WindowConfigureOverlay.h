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
#ifndef WINDOWCONFIGUREOVERLAY_H
#define WINDOWCONFIGUREOVERLAY_H

#include <QWidget>
#include "UIUpdateHelper.h"
#include "WidgetHistogram.h"

namespace Ui
{
class WindowConfigureOverlay;
}

class Layer;
class LayerSurface;
class SurfaceLabel;
class SurfaceOverlayProperty;
class QAbstractButton;
class DialogScreenshotOverlay;

class WindowConfigureOverlay : public QWidget, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit WindowConfigureOverlay(QWidget *parent = 0);
  ~WindowConfigureOverlay();

  void showEvent(QShowEvent *);
  void hideEvent(QHideEvent *);
  void resizeEvent(QResizeEvent* e);

signals:
  void ActiveFrameChanged(int nframe);
  void MaskLoadRequested(const QString& filename);
  void OverlayChanged();

public slots:
  void UpdateGraph(bool bApply = false);
  void UpdateGraphAndApply()
  {
    UpdateGraph(true);
  }
  void UpdateUI();
  void OnCurrentVertexChanged();
  void OnFrameChanged(int nFrame);
  void LoadLabelMask(const QString& fn);

protected slots:
  void OnActiveSurfaceChanged(Layer* layer);
  void OnButtonClicked();
  void OnSliderOpacity( int nVal );
  void OnSpinBoxOpacity( double dVal );
  void OnButtonAdd();
  bool UpdateOverlayProperty( SurfaceOverlayProperty* p );
  void UpdateThresholdChanges();
  void OnHistogramMouseButtonPressed(int button, double value);
  void OnHistogramMarkerChanged();
  void OnSmoothChanged();
  void OnTextThresholdChanged(const QString& strg);
  void OnApply();
  void OnCheckComputeCorrelation(bool bChecked);
  void OnComboCorrelationVolume(int n);
  void OnCheckUsePercentile(bool bChecked);
  void OnCustomColorScale();
  void CheckApply(bool bChecked);
  void OnComboMask(int n);
  void OnCheckInverseMask(bool bChecked);
  void OnSurfaceLabelAdded(SurfaceLabel* label);
  void OnCheckAutoFrameByVertex(bool bChecked);
  void OnCheckUseNonZeroVertices(bool bChecked);
  void OnComboOverlayChanged(int n);
  void OnCycleOverlay();
  void UpdateGeometry();
  void OnCheckFixedAxes(bool bUpdateGraph = true);
  void OnActiveOverlayChanged();
  void OnButtonLoadCustom();
  void OnButtonSaveCustom();
  void OnCheckSetMidToMin(bool bChecked);

private:
  Ui::WindowConfigureOverlay *ui;

  LineMarkers   m_markers;    // custom gradient markers
  LayerSurface* m_layerSurface;
  float*        m_fDataCache;
  double        m_dSavedOffset;
  double        m_rangeOverall[2];
  int           m_nMaxHistCount;
  DialogScreenshotOverlay*  m_dlgScreenshot;
};

#endif // WINDOWCONFIGUREOVERLAY_H
