/**
 * @file  InfoTreeWidget.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/05/29 16:59:22 $
 *    $Revision: 1.8 $
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
#include <QVariantMap>

class QLineEdit;
class QTreeWidgetItem;
class Layer;

class InfoTreeWidget : public QTreeWidget
{
  Q_OBJECT
public:
  InfoTreeWidget(QWidget* parent = 0);
  void contextMenuEvent(QContextMenuEvent *);

signals:
  void RASChangeTriggered(double x, double y, double z);

public slots:
  void UpdateTrackVolumeAnnotation(Layer* layer, const QVariantMap& info);
  void UpdateAll();

protected slots:
  void OnMousePositionChanged();
  void OnCursorPositionChanged();
  void OnItemClicked(QTreeWidgetItem * item, int column);
  void OnEditFinished();
  void OnToggleShowInfo(bool bShow);
  void OnToggleSurfaceCurvature(bool show);
  void OnToggleSurfaceNormal(bool show);

protected:
  void showEvent(QShowEvent *);
  void keyPressEvent(QKeyEvent *event);
  void mousePressEvent(QMouseEvent *event);

private:
  double m_dRAS[3];
  bool  m_bShowSurfaceNormal;
  bool  m_bShowSurfaceCurvature;
  QLineEdit*  m_editor;
  QTreeWidgetItem* m_itemEdited;
};

#endif // INFOTREEWIDGET_H
