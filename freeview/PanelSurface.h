/**
 * @file  PanelSurface.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/04/26 02:38:43 $
 *    $Revision: 1.27 $
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

class PanelSurface : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelSurface(QWidget *parent = 0);
  ~PanelSurface();

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );

protected slots:
  void OnChangeSurfaceType( QAction* act );
  void OnSliderOpacity( int nVal );
  void OnSliderMidPoint( int nVal );
  void OnSliderSlope( int nVal );
  void OnComboCurvature( int nSel );
  void OnLineEditMidPoint( const QString& text );
  void OnLineEditSlope( const QString& text );
  void OnComboOverlay( int nSel );
  void OnComboAnnotation( int nSel );
  void OnComboLabel( int nSel );
  void OnComboVector( int nSel );
  void OnComboSpline(int nSel );
  void OnButtonConfigureOverlay();
  void OnEditPositionOffset();

private:
  Ui::PanelSurface *ui;

  QList<QWidget*>  m_widgetsMidPoint;
  QList<QWidget*>  m_widgetsSlope;
  QList<QWidget*>  m_widgetsVector;
  QList<QWidget*>  m_widgetsVertex;
  QList<QWidget*>  m_widgetsMesh;
  QList<QWidget*>  m_widgetsLabel;
  QList<QWidget*>  m_widgetsSpline;

  WindowConfigureOverlay* m_wndConfigureOverlay;
};

#endif // PANELSURFACE_H
