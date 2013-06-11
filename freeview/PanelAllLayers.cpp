#include "PanelAllLayers.h"
#include "ui_PanelAllLayers.h"
#include "Layer.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include <QDebug>

PanelAllLayers::PanelAllLayers(QWidget *parent) :
  QScrollArea(parent),
  ui(new Ui::PanelAllLayers)
{
  ui->setupUi(this);
  MainWindow* wnd = MainWindow::GetMainWindow();
  QStringList layer_types;
  layer_types << "MRI" << "Surface" << "ROI" << "PointSet" << "CMAT" << "Track";
  foreach (QString type, layer_types)
  {
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerAdded(Layer*)), this, SLOT(OnLayerAdded(Layer*)));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerRemoved(Layer*)), this, SLOT(OnLayerRemoved(Layer*)));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerMoved(Layer*)), this, SLOT(RefreshLayerList()));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerCycled(Layer*)), this, SLOT(RefreshLayerList()));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerVisibilityChanged()), this, SLOT(OnLayerChanged()));
    connect(wnd->GetLayerCollection(type), SIGNAL(ActiveLayerChanged(Layer*)), this, SLOT(OnActiveLayerChanged(Layer*)));
  }
  ui->treeWidgetLayers->setEditTriggers( QTreeWidget::SelectedClicked );
  ui->stackedWidget->hide();

  for (int i = 0; i < ui->stackedWidget->count(); i++)
  {
    PanelLayer* panel = qobject_cast<PanelLayer*>(ui->stackedWidget->widget(i));
    if (panel)
      panel->InitializeLayerTreeWidget(ui->treeWidgetLayers);
  }
}

PanelAllLayers::~PanelAllLayers()
{
  delete ui;
}

void PanelAllLayers::OnActiveLayerChanged(Layer *curlayer)
{
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  if (item)
  {
    QObject* sel = qobject_cast<QObject*>(item->data(0, Qt::UserRole).value<QObject*>());
    if (sel != curlayer)
    {
      for (int i = 0; i < ui->treeWidgetLayers->topLevelItemCount(); i++)
      {
        QTreeWidgetItem* topItem = ui->treeWidgetLayers->topLevelItem(i);
        for (int j = 0; j < topItem->childCount(); j++)
        {
          item = topItem->child(j);
          Layer* layer = qobject_cast<Layer*>(item->data(0, Qt::UserRole).value<QObject*>());
          if (layer == curlayer)
          {
            ui->treeWidgetLayers->setCurrentItem(item);
            return;
          }
        }
      }
    }
  }
}

void PanelAllLayers::OnLayerAdded(Layer *added_layer)
{
  RefreshLayerList(added_layer);
}

void PanelAllLayers::RefreshLayerList(Layer *curLayer)
{
  Layer* layer = NULL;
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  if (item)
    layer = qobject_cast<Layer*>(item->data( 0, Qt::UserRole ).value<QObject*>());
  if (curLayer)
    layer = curLayer;

  ui->treeWidgetLayers->blockSignals(true);
  ui->treeWidgetLayers->clear();
  ui->treeWidgetLayers->blockSignals(false);
  MainWindow* wnd = MainWindow::GetMainWindow();
  AddLayers(wnd->GetLayers("MRI"), "Volumes", layer);
  AddLayers(wnd->GetLayers("Surface"), "Surfaces", layer);
  AddLayers(wnd->GetLayers("ROI"), "ROIs", layer);
  AddLayers(wnd->GetLayers("PointSet"), "Point Sets", layer);
  AddLayers(wnd->GetLayers("CMAT"), "CMAT", layer);
  AddLayers(wnd->GetLayers("Track"), "Track", layer);
}

void PanelAllLayers::AddLayers(QList<Layer *> layers, const QString &cat_name, Layer* curLayer)
{
  ui->treeWidgetLayers->blockSignals(true);
  QTreeWidgetItem* currentItem = NULL;
  if (!layers.isEmpty())
  {
    QTreeWidgetItem* topItem = new QTreeWidgetItem(ui->treeWidgetLayers);
    topItem->setText(0, cat_name);
    topItem->setFlags(topItem->flags() & (~Qt::ItemIsSelectable));
    for (int i = 0; i < layers.size(); i++)
    {
      QTreeWidgetItem* item = new QTreeWidgetItem();
      item->setText(0, layers[i]->GetName());
      item->setData(0, Qt::UserRole, QVariant::fromValue((QObject*)layers[i]));
      item->setCheckState(0, layers[i]->IsVisible() ? Qt::Checked : Qt::Unchecked);
      item->setFlags( item->flags() | Qt::ItemIsEditable );
      topItem->addChild(item);
      if (layers[i] == curLayer)
        currentItem = item;
    }
    topItem->setExpanded(true);
  }
  ui->treeWidgetLayers->blockSignals(false);
  if (currentItem)
    ui->treeWidgetLayers->setCurrentItem(currentItem);
}

void PanelAllLayers::OnLayerRemoved(Layer * removed_layer)
{
  ui->treeWidgetLayers->blockSignals(true);
  for (int i = 0; i < ui->treeWidgetLayers->topLevelItemCount(); i++)
  {
    QTreeWidgetItem* topItem = ui->treeWidgetLayers->topLevelItem(i);
    for (int j = 0; j < topItem->childCount(); j++)
    {
      QTreeWidgetItem* item = topItem->child(j);
      if (item->data(0, Qt::UserRole).value<QObject*>() == removed_layer)
      {
        topItem->removeChild(item);
        delete item;
        if (topItem->childCount() == 0)
        {
          ui->treeWidgetLayers->takeTopLevelItem(ui->treeWidgetLayers->indexOfTopLevelItem(topItem));
          delete topItem;
          topItem = NULL;
        }
        ui->treeWidgetLayers->blockSignals(false);
        if (topItem == NULL)  // removed the last layer in the layer collection, now switch selection to another layer collection if existed.
        {
          if (ui->treeWidgetLayers->topLevelItemCount() > 0)
          {
            topItem = ui->treeWidgetLayers->topLevelItem(0);
            ui->treeWidgetLayers->setCurrentItem(topItem->child(0));
          }
        }
        return;
      }
    }
  }
  ui->treeWidgetLayers->blockSignals(false);
}

PanelLayer* PanelAllLayers::SetCurrentPanel(const QString& layerType)
{
  for (int i = 0; i < ui->stackedWidget->count(); i++)
  {
    PanelLayer* panel = qobject_cast<PanelLayer*>(ui->stackedWidget->widget(i));
    if (panel && panel->GetLayerType() == layerType)
    {
      ui->stackedWidget->setCurrentWidget(panel);
      return panel;
    }
  }
  return NULL;
}

void PanelAllLayers::OnCurrentItemChanged(QTreeWidgetItem *item)
{
  ui->stackedWidget->setVisible(item);
  if (!item)
    return;

  for (int i = 0; i < ui->stackedWidget->count(); i++)
  {
    PanelLayer* panel = qobject_cast<PanelLayer*>(ui->stackedWidget->widget(i));
    if (panel)
      panel->SetCurrentLayer(NULL);
  }

  Layer* layer = qobject_cast<Layer*>(item->data(0, Qt::UserRole).value<QObject*>());
  if (!layer)
    return;

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  QString type = layer->GetPrimaryType();
  LayerCollection* lc = mainwnd->GetLayerCollection(type);
  if (lc)
  {
    lc->SetActiveLayer(layer);
  }

  PanelLayer* panel = SetCurrentPanel(type);
  if (panel)
    panel->SetCurrentLayer(layer);
}

void PanelAllLayers::OnItemChanged(QTreeWidgetItem *item)
{
  Layer* layer = qobject_cast<Layer*>(item->data( 0, Qt::UserRole ).value<QObject*>());
  if ( layer )
  {
    layer->SetName( item->text(0) );
    layer->SetVisible( item->checkState( 0 ) == Qt::Checked );
  }
}

void PanelAllLayers::OnItemSelectionChanged()
{
  // first disconnect all previous layer connections
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  if (item)
  {
    Layer* layer = qobject_cast<Layer*>(item->data(0, Qt::UserRole).value<QObject*>());
    if (layer)
    {
      SetCurrentPanel(layer->GetPrimaryType());
    }
  }
  PanelLayer* panel = qobject_cast<PanelLayer*>(ui->stackedWidget->currentWidget());
  if (panel)
  {
    panel->DisconnectAllLayers();

    QList<Layer*> layers = GetSelectedLayers(panel->GetLayerType());
    if (!layers.isEmpty())
    {
      foreach (Layer* l, layers)
        panel->ConnectLayer(l);
    }
  }
}

void PanelAllLayers::OnLayerChanged()
{
  ui->treeWidgetLayers->blockSignals(true);
  for (int i = 0; i < ui->treeWidgetLayers->topLevelItemCount(); i++)
  {
    QTreeWidgetItem* topItem = ui->treeWidgetLayers->topLevelItem(i);
    for (int j = 0; j < topItem->childCount(); j++)
    {
      QTreeWidgetItem* item = topItem->child(j);
      Layer* layer = qobject_cast<Layer*>(item->data(0, Qt::UserRole).value<QObject*>());
      if (layer)
      {
        item->setText(0, layer->GetName());
        item->setCheckState(0, layer->IsVisible() ?  Qt::Checked : Qt::Unchecked);
      }
    }
  }
  ui->treeWidgetLayers->blockSignals(false);
}

void PanelAllLayers::OnItemDoubleClicked(QTreeWidgetItem *item)
{
  Layer* layer = qobject_cast<Layer*>(item->data( 0, Qt::UserRole ).value<QObject*>());
/*
  QStringList layer_types;
  layer_types << "MRI" << "Surface" << "ROI" << "PointSet" << "CMAT";
  foreach (QString type, layer_types)
  {
    if (layer->IsTypeOf(type))
    {
      emit LayerTypeTriggered(type);
      return;
    }
  }
  */
  if (layer)
  {
    MainWindow* mainwnd = MainWindow::GetMainWindow();
    LayerCollection* lc = mainwnd->GetLayerCollection(GetCurrentLayerType());
    if (lc)
    {
      lc->MoveToTop(layer);
      lc->SetActiveLayer(layer);
    }
  }
}

QString PanelAllLayers::GetCurrentLayerType()
{
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  if (item)
  {
    Layer* layer = qobject_cast<Layer*>(item->data( 0, Qt::UserRole ).value<QObject*>());
    if (layer)
    {
      QStringList layer_types;
      layer_types << "MRI" << "Surface" << "ROI" << "PointSet" << "CMAT" << "Track";
      foreach (QString type, layer_types)
      {
        if (layer->IsTypeOf(type))
        {
          return type;
        }
      }
    }
  }
  return "";
}

QList<Layer*> PanelAllLayers::GetSelectedLayers(const QString &layerType)
{
  QList<Layer*> layers;
  QList<QTreeWidgetItem*> items = ui->treeWidgetLayers->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = qobject_cast<Layer*>(item->data( 0, Qt::UserRole ).value<QObject*>());
    if (layer && layer->IsTypeOf(layerType))
    {
      layers << layer;
    }
  }
  return layers;
}

void PanelAllLayers::UpdateWidgets()
{
  PanelLayer* panel = qobject_cast<PanelLayer*>(ui->stackedWidget->currentWidget());
  if (panel)
    panel->UpdateWidgets();
}
