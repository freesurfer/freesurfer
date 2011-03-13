/**
 * @file  InfoTreeWidget.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:17 $
 *    $Revision: 1.2 $
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

#ifndef INFOTREEWIDGET_H
#define INFOTREEWIDGET_H

#include <QTreeWidget>

class QLineEdit;
class QTreeWidgetItem;

class InfoTreeWidget : public QTreeWidget
{
  Q_OBJECT
public:
  InfoTreeWidget(QWidget* parent = 0);

signals:
  void RASChangeTriggered(double x, double y, double z);

protected slots:
  void OnMousePositionChanged();
  void OnCursorPositionChanged();
  void UpdateAll();
  void OnItemClicked(QTreeWidgetItem * item, int column);
  void OnEditFinished();

protected:
  void showEvent(QShowEvent *);
  void keyPressEvent(QKeyEvent *event);
  void mousePressEvent(QMouseEvent *event);

private:
  double m_dRAS[3];
  QLineEdit*  m_editor;
  QTreeWidgetItem* m_itemEdited;
};

#endif // INFOTREEWIDGET_H
