/**
 * @file  LayerTreeWidget.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2015/10/16 17:31:25 $
 *    $Revision: 1.13 $
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
#include <QItemDelegate>

class Layer;

class MyItemDelegate : public QItemDelegate
{
  Q_OBJECT

public:
  explicit MyItemDelegate (QTreeWidget *parent)
    : QItemDelegate (parent), ParentView (parent) { }
  ~MyItemDelegate() { }

  QRect GetCheckBoxRect(const QModelIndex &index, const QStyleOptionViewItem& option) const;

private:
  QTreeWidget* ParentView ;
};

class LayerTreeWidget : public QTreeWidget
{
  Q_OBJECT
public:
  explicit LayerTreeWidget(QWidget *parent = 0);

  void contextMenuEvent(QContextMenuEvent *e);
  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);

signals:

public slots:
  void ForceUpdate();

  void OnShowAll();
  void OnHideAll();
  void OnLockAll();
  void OnUnlockAll();
  void OnShowAllInfo();
  void OnHideAllInfo();
  void OnSetColorMap();
  void OnEditName();
  void OnSaveVisibleVolumes();
  void SelectAll();
  void DeselectAll();

protected:
  void drawRow ( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const;

  MyItemDelegate* m_itemDelegate;
  QRect         rectCheckbox;
  bool          m_bCheckBoxClicked;
};



#endif // LAYERTREEWIDGET_H
