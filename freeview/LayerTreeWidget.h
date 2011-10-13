/**
 * @file  LayerTreeWidget.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/10/13 21:05:31 $
 *    $Revision: 1.6 $
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
#ifndef LAYERTREEWIDGET_H
#define LAYERTREEWIDGET_H

#include <QTreeWidget>

class LayerTreeWidget : public QTreeWidget
{
  Q_OBJECT
public:
  explicit LayerTreeWidget(QWidget *parent = 0);

  void contextMenuEvent(QContextMenuEvent *e);
  void mousePressEvent(QMouseEvent *event);

signals:

public slots:
  void ForceUpdate();

  void OnShowAll();
  void OnHideAll();
  void OnLockAll();
  void OnUnlockAll();

protected:
  void drawRow ( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const;
};

#endif // LAYERTREEWIDGET_H
