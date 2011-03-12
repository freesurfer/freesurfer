#include "PanelLayer.h"
#include "Layer.h"
#include "LayerProperty.h"
#include "LayerCollection.h"
#include <QTimer>
#include <QApplication>
#include <QTreeWidget>
#include <QAction>
#include <QLineEdit>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QDebug>

PanelLayer::PanelLayer(QWidget *parent) :
    QScrollArea(parent),
    UIUpdateHelper(),
    m_bToUpdate( false ),
    treeWidgetPrivate( NULL )
{
    QTimer* idleTimer = new QTimer(this);
    connect( idleTimer, SIGNAL(timeout()), this, SLOT(OnIdle()), Qt::QueuedConnection);
    idleTimer->start(500);
    QTimer* updateTimer = new QTimer(this);
    connect( updateTimer, SIGNAL(timeout()), this, SLOT(OnUpdate()), Qt::QueuedConnection);
    updateTimer->start(50);
}

void PanelLayer::UpdateWidgets()
{
    m_bToUpdate = true;
}

void PanelLayer::OnIdle()
{
    DoIdle();
}

void PanelLayer::OnUpdate()
{
    if ( m_bToUpdate )
    {
        DoUpdateWidgets();
        m_bToUpdate = false;
    }
}

void PanelLayer::InitializeLayerList( QTreeWidget* treeWidget, LayerCollection* cl )
{
    allWidgets = this->widget()->findChildren<QWidget*>();
    allActions = this->findChildren<QAction*>();
    treeWidgetPrivate = treeWidget;
    layerCollectionPrivate = cl;
//    treeWidget->setEditTriggers( QTreeWidget::SelectedClicked );
    connect( cl, SIGNAL(LayerAdded(Layer*)), this, SLOT(OnLayerAdded(Layer*)), Qt::UniqueConnection );
    connect( cl, SIGNAL(LayerRemoved(Layer*)), this, SLOT(OnLayerRemoved(Layer*)), Qt::UniqueConnection );
    connect( cl, SIGNAL(LayerMoved(Layer*)), this, SLOT(OnLayerMoved(Layer*)), Qt::UniqueConnection );
    connect( cl, SIGNAL(ActiveLayerChanged(Layer*)), this, SLOT(OnActiveLayerChanged(Layer*)), Qt::UniqueConnection );
    connect( cl, SIGNAL(LayerCycled(Layer*)), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
    connect( cl, SIGNAL(LayerVisibilityChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
    connect( treeWidget, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)),
             this, SLOT(OnCurrentItemChanged(QTreeWidgetItem*)) );
    connect( treeWidget, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)),
             this, SLOT(UpdateWidgets()));
    connect( treeWidget, SIGNAL(itemChanged(QTreeWidgetItem*,int)), this, SLOT(OnItemChanged(QTreeWidgetItem*)) );
    connect( treeWidget, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(OnItemDoubleClicked(QModelIndex)));

    UpdateWidgets();
}

void PanelLayer::ConnectLayer( Layer* layer )
{
  connect( layer, SIGNAL(Locked(bool)), treeWidgetPrivate, SLOT(repaint()) );
}

void PanelLayer::OnLayerAdded( Layer* layer )
{
    QTreeWidget* t = treeWidgetPrivate;
    for ( int i = 0; i < t->topLevelItemCount(); i++ )
    {
        QTreeWidgetItem* item = t->topLevelItem( i );
        if ( item->data( 0, Qt::UserRole ).value<QObject*>() == layer )
        {
            t->setCurrentItem( item );
            return;
        }
    }
    t->blockSignals( true );
    QTreeWidgetItem* item = new QTreeWidgetItem();
    item->setData( 0, Qt::UserRole, QVariant::fromValue((QObject*)layer) );
    item->setText( 0, layer->GetName() );
    item->setCheckState( 0, layer->IsVisible() ? Qt::Checked : Qt::Unchecked );
    item->setFlags( item->flags() | Qt::ItemIsEditable );
    t->insertTopLevelItem( 0, item );
    t->blockSignals( false );
    t->setCurrentItem( item );
}

void PanelLayer::OnLayerRemoved( Layer *layer )
{
    QTreeWidget* t = treeWidgetPrivate;
    for ( int i = 0; i < t->topLevelItemCount(); i++ )
    {
        QTreeWidgetItem* item = t->topLevelItem( i );
        if ( item->data( 0, Qt::UserRole ).value<QObject*>() == layer )
        {
            t->takeTopLevelItem( i );
            return;
        }
    }
}

void PanelLayer::OnLayerMoved(Layer *layer_in)
{
  QTreeWidget* t = treeWidgetPrivate;
  t->blockSignals( true );
  for ( int i = 0; i < t->topLevelItemCount(); i++ )
  {
      QTreeWidgetItem* item = t->topLevelItem( i );
      Layer* layer = layerCollectionPrivate->GetLayer( i );
      item->setData( 0, Qt::UserRole, QVariant::fromValue((QObject*)layer) );
      item->setText( 0, layer->GetName() );
      item->setCheckState( 0, layer->IsVisible() ? Qt::Checked : Qt::Unchecked );
      if ( layer == layer_in )
      {
        t->setCurrentItem( item );
      }
  }
  t->blockSignals( false );
}

void PanelLayer::OnActiveLayerChanged(Layer *layer_in)
{
    QTreeWidget* t = treeWidgetPrivate;
    for ( int i = 0; i < t->topLevelItemCount(); i++ )
    {
        QTreeWidgetItem* item = t->topLevelItem( i );
        if ( item->data( 0, Qt::UserRole ).value<QObject*>() == layer_in )
        {
            t->setCurrentItem( item );
            break;
        }
    }
}

void PanelLayer::OnCurrentItemChanged(QTreeWidgetItem *item)
{
    if (item)
    {
        Layer* layer = qobject_cast<Layer*>(item->data( 0, Qt::UserRole ).value<QObject*>());
        if (layer)
        {
            layerCollectionPrivate->SetActiveLayer( layer );
            // first disconnect all previous layer connections
            for ( int i = 0; i < layerCollectionPrivate->GetNumberOfLayers(); i++ )
            {
              for ( int j = 0; j < allWidgets.size(); j++ )
              {
                layerCollectionPrivate->GetLayer( i )->disconnect( this );
                layerCollectionPrivate->GetLayer( i )->GetProperty()->disconnect( this );
                allWidgets[j]->disconnect( layerCollectionPrivate->GetLayer( i ) );
                allWidgets[j]->disconnect( layerCollectionPrivate->GetLayer( i )->GetProperty() );
              }
            }
            ConnectLayer( layer );
        }
    }
}

void PanelLayer::OnItemChanged( QTreeWidgetItem* item )
{
  Layer* layer = qobject_cast<Layer*>(item->data( 0, Qt::UserRole ).value<QObject*>());
  if ( layer )
  {
    layer->SetName( item->text(0) );
    layer->SetVisible( item->checkState( 0 ) == Qt::Checked );
  }
}

void PanelLayer::OnItemDoubleClicked(const QModelIndex &index)
{
    Layer* layer = qobject_cast<Layer*>(index.data(Qt::UserRole).value<QObject*>());
    if (layer)
    {
        layerCollectionPrivate->MoveToTop(layer);
        layerCollectionPrivate->SetActiveLayer(layer);
    }
}

void PanelLayer::BlockAllSignals( bool bBlock )
{
    for ( int i = 0; i < allWidgets.size(); i++ )
        allWidgets[i]->blockSignals( bBlock );
    for ( int i = 0; i < allActions.size(); i++ )
        allActions[i]->blockSignals( bBlock );
}

