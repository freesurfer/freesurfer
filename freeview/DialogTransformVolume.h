/**
 * @file  DialogTransformVolume.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
 *    $Revision: 1.10 $
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
#ifndef DIALOGTRANSFORMVOLUME_H
#define DIALOGTRANSFORMVOLUME_H

#include <QDialog>
#include "UIUpdateHelper.h"

namespace Ui
{
class DialogTransformVolume;
}

class QCheckBox;
class QComboBox;
class QLineEdit;
class QScrollBar;

class DialogTransformVolume : public QDialog, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit DialogTransformVolume(QWidget *parent = 0);
  ~DialogTransformVolume();

  bool GetRotation( int nID_in, int& plane_out, double& angle_out );
  void UpdateUI( int scope = 2 );

protected slots:
  void OnApply();
  void OnRestore();
  void OnSaveReg();

  void OnScrollBarTranslateX(int nVal);
  void OnScrollBarTranslateY(int nVal);
  void OnScrollBarTranslateZ(int nVal);
  void OnLineEditTranslateX(const QString& text);
  void OnLineEditTranslateY(const QString& text);
  void OnLineEditTranslateZ(const QString& text);
  void OnScrollBarScaleX(int nVal);
  void OnScrollBarScaleY(int nVal);
  void OnScrollBarScaleZ(int nVal);
  void OnLineEditScaleX(const QString& text);
  void OnLineEditScaleY(const QString& text);
  void OnLineEditScaleZ(const QString& text);

  void OnActiveLayerChanged();

private:
  void DoRotate();
  void RespondTextTranslate   ( int n );
  void RespondScrollTranslate ( int n );
  void RespondTextScale   ( int n );
  void RespondScrollScale ( int n );

  Ui::DialogTransformVolume *ui;

  QCheckBox*   m_checkRotate[3];
  QComboBox*   m_comboRotate[3];
  QLineEdit*   m_textAngle[3];
  QScrollBar*    m_scrollTranslate[3];
  QLineEdit*     m_textTranslate[3];
  QScrollBar*    m_scrollScale[3];
  QLineEdit*     m_textScale[3];
};

#endif // DIALOGTRANSFORMVOLUME_H
