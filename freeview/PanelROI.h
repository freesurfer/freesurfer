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
#ifndef PANELROI_H
#define PANELROI_H

#include "PanelLayer.h"
#include <QList>

namespace Ui
{
class PanelROI;
}

class PanelROI : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelROI(QWidget *parent = 0);
  ~PanelROI();

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );

protected slots:
  void OnSliderOpacity( int val );
  void OnEditThreshold(const QString& text);
  void OnEditHeatscaleMin(const QString& text);
  void OnEditHeatscaleMax(const QString& text);
  void OnComboMappedSurface(int nIndex);
  void OnButtonDilate();
  void OnButtonErode();
  void OnButtonOpen();
  void OnButtonClose();
  void OnButtonResample();

private:
  Ui::PanelROI *ui;
  QList<QWidget*>   m_listWidgetsHeatscale;
};

#endif // PANELROI_H
