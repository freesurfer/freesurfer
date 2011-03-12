#ifndef WINDOWCONFIGUREOVERLAY_H
#define WINDOWCONFIGUREOVERLAY_H

#include <QWidget>
#include "UIUpdateHelper.h"
#include "WidgetHistogram.h"

namespace Ui {
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

private:
    Ui::WindowConfigureOverlay *ui;

    LineMarkers   m_markers;    // custom gradient markers
    LayerSurface* m_layerSurface;
};

#endif // WINDOWCONFIGUREOVERLAY_H
