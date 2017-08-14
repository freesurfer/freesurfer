/**
 * @file  LayerTreeWidget.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2016/05/31 18:30:40 $
 *    $Revision: 1.15 $
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
class QDropEvent;

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
  void ToReorderLayers(const QList<Layer*>& newlist);

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
  void selectAll()
  {
    SelectAll();
  }
  void DeselectAll();

protected:
  bool event(QEvent* e);
  void drawRow ( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const;
  virtual void dropEvent(QDropEvent * event);

  MyItemDelegate* m_itemDelegate;
  QRect         rectCheckbox;
  bool          m_bCheckBoxClicked;
};



#endif // LAYERTREEWIDGET_H
