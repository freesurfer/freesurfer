/**
 * @file  PanelTrack.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/04/03 21:07:22 $
 *    $Revision: 1.4.2.1 $
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
#ifndef PANELTRACK_H
#define PANELTRACK_H

#include "PanelLayer.h"
#include <QList>

namespace Ui
{
class PanelTrack;
}

class PanelTrack : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelTrack(QWidget *parent = 0);
  ~PanelTrack();

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );

private:
  Ui::PanelTrack *ui;

  QList<QWidget*> m_widgetlistDirectionalColor;
  QList<QWidget*> m_widgetlistSolidColor;
};

#endif // PANELTRACK_H
