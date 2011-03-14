/**
 * @file  PanelPointSet.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
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
#ifndef PANELPOINTSET_H
#define PANELPOINTSET_H

#include "PanelLayer.h"

namespace Ui
{
class PanelPointSet;
}

class PanelPointSet : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelPointSet(QWidget *parent = 0);
  ~PanelPointSet();

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );
  void LoadScalarValues();

protected slots:
  void OnSliderOpacity( int nVal );
  void OnSliderMin(int nVal);
  void OnSliderMid(int nVal);
  void OnSliderMax(int nVal);
  void OnSliderOffset(int nVal);
  void OnLineEditMin(const QString& text);
  void OnLineEditMid(const QString& text);
  void OnLineEditMax(const QString& text);
  void OnLineEditOffset(const QString& text);
  void OnLineEditRadius(const QString& text);
  void OnLineEditSplineRadius(const QString& text);
  void OnComboScalarMap(int nSel);

private:
  Ui::PanelPointSet *ui;
  QList<QWidget*> m_widgetlistSolidColor;
  QList<QWidget*> m_widgetlistHeatScale;
  QList<QWidget*> m_widgetlistSpline;
};

#endif // PANELPOINTSET_H
