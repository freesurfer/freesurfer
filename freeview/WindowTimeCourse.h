/**
 * @file  WindowTimeCourse.h
 * @brief Tool window to display time course data
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2014/04/09 20:56:04 $
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

#ifndef WINDOWTIMECOURSE_H
#define WINDOWTIMECOURSE_H

#include <QWidget>

namespace Ui {
class WindowTimeCourse;
}

class LayerMRI;
class LayerSurface;
class SurfaceOverlay;

class WindowTimeCourse : public QWidget
{
  Q_OBJECT

public:
  explicit WindowTimeCourse(QWidget *parent = 0);
  ~WindowTimeCourse();

  void showEvent(QShowEvent* e);

public slots:
  void UpdateData(bool bForce = false);
  void OnFrameChanged(int n);
  void SetCurrentFrame(int n);
  void OnLayerCorrelationSurfaceChanged();
  void OnLineEditScaleReturnPressed();
  void OnCheckAutoScale(bool bChecked);
  void OnCheckMaxScale(bool bChecked);
  void OnComboSecondVolume(int nSel);
  void UpdateScaleInfo();
  void UpdateUI();
  void OnCheckShowFrameNumber(bool);

signals:
  void FrameChanged(int frame);

private:

  Ui::WindowTimeCourse *ui;
  LayerMRI*     lastMRI;
  LayerSurface* lastSurface;
  SurfaceOverlay* lastOverlay;
};

#endif // WINDOWTIMECOURSE_H
