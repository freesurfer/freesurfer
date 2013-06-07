/**
 * @file  PanelLayer.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2013/06/07 02:20:32 $
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
#ifndef PANELLAYER_H
#define PANELLAYER_H

#include "UIUpdateHelper.h"
#include <QScrollArea>
#include <QList>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include "LayerCollection.h"

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
  explicit PanelLayer(const QString& layerType, QWidget *parent = 0);

  void InitializeLayerTreeWidget( QTreeWidget* treeWidget);
  void SetCurrentLayer(Layer* layer);
  void DisconnectAllLayers();
  virtual void ConnectLayer( Layer* layer );

  QString GetLayerType()
  {
    return m_layerType;
  }

signals:

public slots:
  void UpdateWidgets();

protected:
  virtual void DoIdle() = 0;
  virtual void DoUpdateWidgets() = 0;
  void BlockAllSignals( bool bBlock );
  template<typename T> inline T GetCurrentLayer();
  template<typename T> inline QList<T> GetSelectedLayers();

protected slots:
  void OnIdle();
  void OnUpdate();
//  void OnLayerMoved( Layer* layer );
//  void OnActiveLayerChanged( Layer* layer );

  void ResetLayerCollection()
  {
    m_layerCollection = NULL;
  }

protected:
  QList<QWidget*>     allWidgets;
  QList<QAction*>     allActions;
  QTreeWidget*        treeWidgetLayers;
  LayerCollection*    m_layerCollection;

private:
  bool m_bToUpdate;
  Layer*              m_currentLayer;
  QString             m_layerType;
};

template<typename T> T PanelLayer::GetCurrentLayer()
{
  if ( !treeWidgetLayers || !m_layerCollection)
  {
    return NULL;
  }

  if (m_layerCollection->Contains(m_currentLayer))
    return qobject_cast<T>(m_currentLayer);
  else
    return NULL;
}

template<typename T> QList<T> PanelLayer::GetSelectedLayers()
{
  QList<T> list;
  if ( !treeWidgetLayers )
  {
    return list;
  }

  QList<QTreeWidgetItem*> items = treeWidgetLayers->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    T t = qobject_cast<T>( item->data(0, Qt::UserRole).template value<QObject*>() );
    if (t)
      list << t;
  }
  return list;
}

#endif // PANELLAYER_H
