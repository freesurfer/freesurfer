/**
 * @file  PanelLayer.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/03/13 21:32:06 $
 *    $Revision: 1.7 $
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
#ifndef PANELLAYER_H
#define PANELLAYER_H

#include "UIUpdateHelper.h"
#include <QScrollArea>
#include <QList>
#include <QTreeWidget>
#include <QTreeWidgetItem>

class QLineEdit;
class QSpinBox;
class QDoubleSpinBox;
class QModelIndex;
class Layer;
class LayerCollection;

class PanelLayer : public QScrollArea, public UIUpdateHelper
{
  Q_OBJECT
public:
  explicit PanelLayer(QWidget *parent = 0);

signals:

public slots:
  void UpdateWidgets();

protected:
  virtual void DoIdle() = 0;
  virtual void DoUpdateWidgets() = 0;
  virtual void ConnectLayer( Layer* layer );
  void InitializeLayerList( QTreeWidget* treeWidget, LayerCollection* lc );
  void BlockAllSignals( bool bBlock );
  template<typename T> inline T GetCurrentLayer();
  template<typename T> inline QList<T> GetSelectedLayers();

protected slots:
  void OnIdle();
  void OnUpdate();
  void OnLayerAdded( Layer* layer );
  void OnLayerRemoved( Layer* layer );
  void OnLayerMoved( Layer* layer );
  void OnActiveLayerChanged( Layer* layer );
  void OnItemChanged( QTreeWidgetItem* item );
  void OnCurrentItemChanged( QTreeWidgetItem* item );
  void OnItemSelectionChanged();
  void OnItemDoubleClicked(const QModelIndex& index);
  void OnLayerNameChanged();

protected:
  QList<QWidget*>     allWidgets;
  QList<QAction*>     allActions;

private:
  bool m_bToUpdate;
  QTreeWidget*        treeWidgetPrivate;
  LayerCollection*    layerCollectionPrivate;
};

template<typename T> T PanelLayer::GetCurrentLayer()
{
  if ( !treeWidgetPrivate )
  {
    return NULL;
  }

  QTreeWidgetItem* item = treeWidgetPrivate->currentItem();
  if ( item )
  {
    return qobject_cast<T>( item->data(0, Qt::UserRole).template value<QObject*>() );
  }
  else
  {
    return NULL;
  }
}

template<typename T> QList<T> PanelLayer::GetSelectedLayers()
{
  QList<T> list;
  if ( !treeWidgetPrivate )
  {
    return list;
  }

  QList<QTreeWidgetItem*> items = treeWidgetPrivate->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    list << qobject_cast<T>( item->data(0, Qt::UserRole).template value<QObject*>() );
  }
  return list;
}
#endif // PANELLAYER_H
