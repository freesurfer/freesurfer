/**
 * @file  WindowConfigureOverlay.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/11/01 19:21:06 $
 *    $Revision: 1.8 $
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
class SurfaceOverlayProperty;
class QAbstractButton;

class WindowConfigureOverlay : public QWidget, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit WindowConfigureOverlay(QWidget *parent = 0);
  ~WindowConfigureOverlay();

  virtual void showEvent(QShowEvent *);

signals:
  void ActiveFrameChanged();

public slots:
  void UpdateGraph();
  void UpdateUI();

protected slots:
  void OnActiveSurfaceChanged(Layer* layer);
  void OnClicked( QAbstractButton* btn );
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
  void OnFrameChanged(int nFrame);

private:
  Ui::WindowConfigureOverlay *ui;

  LineMarkers   m_markers;    // custom gradient markers
  LayerSurface* m_layerSurface;
  float*        m_fDataCache;
};

#endif // WINDOWCONFIGUREOVERLAY_H
