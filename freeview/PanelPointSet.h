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
#ifndef PANELPOINTSET_H
#define PANELPOINTSET_H

#include "PanelLayer.h"
#include <QMap>

namespace Ui
{
class PanelPointSet;
}

class QLabel;
class QTreeWidgetItem;

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
  void UpdatePointInfo();

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
  void OnSpinBoxGoToPoint(int val);
  void OnButtonGoToPoint();
  void OnButtonCommentAdd();
  void OnButtonStatAdd();
  void OnButtonStatDelete();
  void OnCommentLabelClicked(const QString& link);
  void ScrollCommentsToBottom();
  void OnStatItemChanged(QTreeWidgetItem* item, int col);
  void OnCurrentStatItemChanged(QTreeWidgetItem* cur, QTreeWidgetItem* old);
  void SetCurrentPoint(int nIndex);
  void OnTextOverallQualityChanged();
  void OnSpinBoxOverallScore(int);
  void OnSpinBoxSecondQA(int);

private:
  QLabel* MakeCommentItem(const QVariantMap& map, QLabel* label_in = NULL);
  QTreeWidgetItem* AddStatItem(const QString& name, double value);

  Ui::PanelPointSet *ui;
  QList<QWidget*> m_widgetlistSolidColor;
  QList<QWidget*> m_widgetlistHeatScale;
  QList<QWidget*> m_widgetlistSpline;

  QString     m_self;
  QMap<QObject*, int> m_mapCurrentPoint;
};

#endif // PANELPOINTSET_H
