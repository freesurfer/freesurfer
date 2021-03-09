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
#ifndef DIALOGTRANSFORMVOLUME_H
#define DIALOGTRANSFORMVOLUME_H

#include <QDialog>
#include "UIUpdateHelper.h"
#include <QList>
#include <QIcon>

namespace Ui
{
class DialogTransformVolume;
}

class QCheckBox;
class QComboBox;
class QLineEdit;
class QScrollBar;
class QtColorPicker;
class QPushButton;
class QSlider;

class DialogTransformVolume : public QDialog, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit DialogTransformVolume(QWidget *parent = 0);
  ~DialogTransformVolume();

  bool GetRotation( int nID_in, int& plane_out, double& angle_out );
  void UpdateUI( int scope = 2 );

  void closeEvent(QCloseEvent * e);
  void showEvent(QShowEvent * e);

signals:
  void CurrentLandmarkChanged(int n);

protected slots:
  void OnApply();
  void OnRestore();
  void OnSaveReg();

  void OnSliderRotateX(int nVal);
  void OnSliderRotateY(int nVal);
  void OnSliderRotateZ(int nVal);
  void OnLineEditRotateX(const QString& text);
  void OnLineEditRotateY(const QString& text);
  void OnLineEditRotateZ(const QString& text);
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

  void OnCheckBoxFlip();

  void OnSampleMethodChanged();

  void OnActiveLayerChanged();

  void OnRadioButtonLandmark(bool bChecked);

  void OnButtonLandmarkPick();

  void UpdateLandmarkColors();

  void OnButtonCenterToCursor();

  void OnCheckBoxApplyToAll(bool bAll);

private:
  void DoRotate();
  void RespondTextTranslate   ( int n );
  void RespondScrollTranslate ( int n );
  void RespondTextScale   ( int n );
  void RespondScrollScale ( int n );
  void RespondTextRotate  (int n);
  void RespondSliderRotate  (int n);

  QIcon MakeIcon(const QColor& color, int size);

  Ui::DialogTransformVolume *ui;

  QCheckBox*   m_checkRotate[3];
  QComboBox*   m_comboRotate[3];
  QScrollBar*     m_sliderRotate[3];
  QLineEdit*   m_textAngle[3];
  QScrollBar*    m_scrollTranslate[3];
  QLineEdit*     m_textTranslate[3];
  QScrollBar*    m_scrollScale[3];
  QLineEdit*     m_textScale[3];
  QList<QtColorPicker*> m_colorPickerLandmark;
  QList<QPushButton*>   m_btnPickLandmark;
  QList<QComboBox*> m_comboLandmark;
  double        m_dIncrementRotate;
  double        m_dIncrementScale;
  double        m_dIncrementTranslate;
};

#endif // DIALOGTRANSFORMVOLUME_H
