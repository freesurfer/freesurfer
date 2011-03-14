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

protected slots:
    void OnIdle();
    void OnUpdate();
    void OnLayerAdded( Layer* layer );
    void OnLayerRemoved( Layer* layer );
    void OnLayerMoved( Layer* layer );
    void OnActiveLayerChanged( Layer* layer );
    void OnItemChanged( QTreeWidgetItem* item );
    void OnCurrentItemChanged( QTreeWidgetItem* item );
    void OnItemDoubleClicked(const QModelIndex& index);

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
    return NULL;

  QTreeWidgetItem* item = treeWidgetPrivate->currentItem();
  if ( item )
    return qobject_cast<T>( item->data(0, Qt::UserRole).template value<QObject*>() );
  else
    return NULL;
}
#endif // PANELLAYER_H
