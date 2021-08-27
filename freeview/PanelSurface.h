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
#ifndef PANELSURFACE_H
#define PANELSURFACE_H

#include "PanelLayer.h"
#include <QList>

namespace Ui
{
class PanelSurface;
}

class QAction;
class WindowConfigureOverlay;
class SurfaceLabel;
class SurfaceSpline;
class QToolButton;
class QActionGroup;
class DialogCustomFill;
class DialogSurfaceLabelOperations;
class WindowEditAnnotation;

class PanelSurface : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelSurface(QWidget *parent = 0);
  ~PanelSurface();

  bool eventFilter(QObject *watched, QEvent *event);

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );
  virtual void DisconnectAllLayers();

protected slots:
  void OnChangeSurfaceType( QAction* act );
  void OnSliderOpacity( int nVal );
  void OnSliderMidPoint( int nVal );
  void OnSliderSlope( int nVal );
  void OnComboCurvature( int nSel );
  void OnLineEditMidPoint( const QString& text );
  void OnLineEditSlope( const QString& text );
  void OnLineEditLabelThreshold(const QString& text);
  void OnLineEditLabelHeatscaleMin(const QString& text);
  void OnLineEditLabelHeatscaleMax(const QString& text);
  void OnComboLabelColorCode(int nSel);
  void OnComboOverlay( int nSel );
  void OnComboAnnotation( int nSel );
//  void OnComboLabel( int nSel );
  void OnComboVector( int nSel );
  void OnComboSpline(int nSel );
  void OnComboColor(int nSel);
  void OnButtonConfigureOverlay();
  void OnButtonEditAnnotation();
  void OnButtonRemoveOverlay();
  void OnEditPositionOffset();
  void OnLabelItemChanged(QTreeWidgetItem *item);
  void OnLabelItemDoubleClicked(QTreeWidgetItem* item);
  void OnCurrentLabelItemChanged(QTreeWidgetItem *item);
  void OnButtonLoadLabel();
  void OnButtonDeleteLabel();
  void OnButtonNewLabel();
  void OnButtonSaveLabel();
  void OnSaveLabelAs();
  void OnToggleOverlay(bool bShow);
  void OnToggleAnnotation(bool bShow);
  void OnColorPickerLabelColor(const QColor& color);
  void OnCheckBoxLabelOutline(bool outline);
  void UpdateLabelWidgets(bool block_signals = true);
  void UpdateSplineWidgets();
  void OnLockLayer(bool b);
  void OnButtonLoadSpline();
  void OnButtonDeleteSpline();
  void OnColorPickerSplineColor(const QColor&);
  void OnCheckBoxSplineProjection(bool);
  void OnSplineItemChanged(QTreeWidgetItem *item);
  void OnCurrentSplineItemChanged(QTreeWidgetItem *item);
  void OnButtonEditCut(bool b);
  void OnButtonEditPath(bool b);
  void OnButtonCutLine();
  void OnButtonCutClosedLine();
  void OnButtonClearCuts();
  void OnButtonFillUncutArea();
  void OnButtonUndoCut();
  void OnLabelResample();
  void OnLabelMaskOverlay();
  void OnLabelMoreOps();
  void OnLabelOperation(const QVariantMap& op);
  void OnSpinBoxZOrder(int nOrder);
  void OnButtonMakePath();
  void OnButtonMakeClosedPath();
  void OnButtonDeletePath();
  void OnButtonCustomFillPath();
  void OnCustomFillTriggered(const QVariantMap& options);
  void OnButtonClearMarks();
  void OnLineEditLabelOpacity(const QString& text);
  void OnButtonLabelUp();
  void OnButtonLabelDown();
  void SetOverlayFrame(int nFrame);
  void OnButtonSaveAnnotation();
  void OnCycleAnnotation();

private:
  QList<SurfaceLabel*> GetSelectedLabels();
  QList<SurfaceSpline*> GetSelectedSplines();

  Ui::PanelSurface *ui;

  QList<QWidget*>  m_widgetsMidPoint;
  QList<QWidget*>  m_widgetsSlope;
  QList<QWidget*>  m_widgetsVector;
  QList<QWidget*>  m_widgetsVertex;
  QList<QWidget*>  m_widgetsMesh;
  QList<QWidget*>  m_widgetsLabel;
  QList<QWidget*>  m_widgetsSpline;
  QList<QWidget*>  m_widgetsOverlay;
  QList<QWidget*>  m_widgetsAnnotation;
  QToolButton*     m_toolButtonSurface;
  QActionGroup*    m_actGroupSurface;

  WindowConfigureOverlay* m_wndConfigureOverlay;
  DialogCustomFill*     m_dlgCustomFill;
  DialogSurfaceLabelOperations* m_dlgLabelOps;
  WindowEditAnnotation*   m_wndEditAnnotation;
};

#endif // PANELSURFACE_H
