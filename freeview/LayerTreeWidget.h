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
#ifndef LAYERTREEWIDGET_H
#define LAYERTREEWIDGET_H

#include <QTreeWidget>
#include <QItemDelegate>
#include <QList>

class Layer;
class LayerMRI;
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

  QList<LayerMRI*> GetLinkedVolumes()
  {
      return m_linkedVolumes;
  }

signals:
  void ToReorderLayers(const QList<Layer*>& newlist);

public slots:
  void ForceUpdate();

  void OnShowAll();
  void OnHideAll();
  void OnLockAll();
  void OnUnlockAll();
  void OnLockOthers();
  void OnUnlockOthers();
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
  void SetSelectedLayers(const QList<int>& layer_ids);
  void OnLinkVolumes();
  void OnUnlinkVolumes();
  void LinkVolume(LayerMRI* vol);

protected:
  bool event(QEvent* e);
  void drawRow ( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const;
  virtual void dropEvent(QDropEvent * event);

  MyItemDelegate* m_itemDelegate;
  QRect         rectCheckbox;
  bool          m_bCheckBoxClicked;
  QList<LayerMRI*>  m_linkedVolumes;
};



#endif // LAYERTREEWIDGET_H
