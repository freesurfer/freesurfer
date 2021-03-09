/**
 * @brief Tool window to display time course data
 *
 */
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

#ifndef WINDOWTIMECOURSE_H
#define WINDOWTIMECOURSE_H

#include <QWidget>
#include <QLabel>

namespace Ui {
class WindowTimeCourse;
}

class LayerMRI;
class LayerSurface;
class SurfaceOverlay;
class FlowLayout;
struct TimeCourseData;

class ClickableLabel : public QLabel
{
    Q_OBJECT

public:
    explicit ClickableLabel(QWidget* parent = Q_NULLPTR, Qt::WindowFlags f = Qt::WindowFlags());
    ~ClickableLabel();

signals:
    void clicked();

protected:
    void mousePressEvent(QMouseEvent* event);

};

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
  void OnComboSecondPlot(int nSel);
  void UpdateScaleInfo();
  void UpdateUI();
  void OnCheckShowFrameNumber(bool);
  void UpdateAll()
  {
    UpdateUI();
    UpdateData();
  }

  void Clear();
  void OnCheckBoxShowData(bool bShow);
  void OnLegendLabelClicked();

signals:
  void OverlayFrameChanged(int frame);

private:
  QWidget* MakeLegendWidget(QObject* obj, const TimeCourseData& td);

  Ui::WindowTimeCourse *ui;
  LayerMRI*     lastMRI;
  LayerSurface* lastSurface;
  SurfaceOverlay* lastOverlay;
  FlowLayout*   layoutLegend;
};

#endif // WINDOWTIMECOURSE_H
