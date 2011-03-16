/**
 * @file  ToolWindowEdit.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/16 22:07:51 $
 *    $Revision: 1.18 $
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
#ifndef TOOLWINDOWEDIT_H
#define TOOLWINDOWEDIT_H

#include "UIUpdateHelper.h"
#include <QWidget>
#include <QList>

namespace Ui
{
class ToolWindowEdit;
}

class ToolWindowEdit : public QWidget, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit ToolWindowEdit(QWidget *parent = 0);
  ~ToolWindowEdit();

public slots:
  void UpdateWidgets();

protected slots:
  void OnIdle();
  void OnEditMode( QAction* act );
  void OnComboReference(int sel);
  void OnLineEditContourValue(const QString& strg);
  void OnLineEditSmoothSD(const QString& strg);
  void OnDrawRangeChanged(const QString& strg);
  void OnExcludeRangeChanged(const QString& strg);
  void OnReplaceLabel();

protected:
  virtual void showEvent(QShowEvent *);

private:
  Ui::ToolWindowEdit *ui;

  bool m_bToUpdateWidgets;
  QList<QWidget*>  m_widgetsBrushSize;
  QList<QWidget*>  m_widgetsReference;
  QList<QWidget*>  m_widgetsTolerance;
  QList<QWidget*>  m_widgetsConstrain;
  QList<QWidget*>  m_widgetsSmooth;
  QList<QWidget*>  m_widgetsContour;
};

#endif // TOOLWINDOWEDIT_H
