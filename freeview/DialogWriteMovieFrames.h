/**
 * @file  DialogWriteMovieFrames.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2013/11/05 20:25:28 $
 *    $Revision: 1.9 $
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
#ifndef DIALOGWRITEMOVIEFRAMES_H
#define DIALOGWRITEMOVIEFRAMES_H

#include <QDialog>
#include <QTimer>

namespace Ui
{
class DialogWriteMovieFrames;
}

class RenderView;

class DialogWriteMovieFrames : public QDialog
{
  Q_OBJECT

public:
  explicit DialogWriteMovieFrames(QWidget *parent = 0);
  ~DialogWriteMovieFrames();

  void showEvent(QShowEvent *e);
  void closeEvent(QCloseEvent *e);

public slots:
  void UpdateUI(bool bUpdateSliceNumber = true);

signals:
  void Started();
  void Stopped();
  void Progress(int n);

protected slots:
  void OnWrite();
  void OnAbort();
  void OnOpen();

  void OnTimeOut();

private:
  Ui::DialogWriteMovieFrames *ui;
  RenderView* m_view;
  QTimer      m_timer;
  double      m_dAngleStepSize;
  int         m_nStartSlice;
  int         m_nSliceStepSize;
  QString     m_strOutputDir;
  int         m_nStepCount;
  int         m_nTotalSteps;
  bool        m_b3D;
};

#endif // DIALOGWRITEMOVIEFRAMES_H
