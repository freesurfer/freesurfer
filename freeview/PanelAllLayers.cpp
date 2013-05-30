#include "PanelAllLayers.h"
#include "ui_PanelAllLayers.h"
#include "Layer.h"
#include "MainWindow.h"
#include "LayerCollection.h"

PanelAllLayers::PanelAllLayers(QWidget *parent) :
  QScrollArea(parent),
  ui(new Ui::PanelAllLayers)
{
  ui->setupUi(this);
  MainWindow* wnd = MainWindow::GetMainWindow();
  QStringList layer_types;
  layer_types << "MRI" << "Surface" << "ROI" << "PointSet" << "CMAT";
  foreach (QString type, layer_types)
  {
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerAdded(Layer*)), this, SLOT(OnLayerAdded(Layer*)));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerRemoved(Layer*)), this, SLOT(OnLayerRemoved(Layer*)));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerMoved(Layer*)), this, SLOT(RefreshLayerList()));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerCycled(Layer*)), this, SLOT(RefreshLayerList()));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerVisibilityChanged()), this, SLOT(OnLayerChanged()));
  }
  connect(this, SIGNAL(LayerTypeTriggered(QString)), wnd, SLOT(SetActivePanel(QString)));
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
    Layer* sel = qobject_cast<Layer*>(item->data(0, Qt::UserRole).value<QObject*>());
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
  ui->treeWidgetLayers->clear();
  MainWindow* wnd = MainWindow::GetMainWindow();
  AddLayers(wnd->GetLayers("MRI"), "Volumes");
  AddLayers(wnd->GetLayers("Surface"), "Surfaces");
  AddLayers(wnd->GetLayers("ROI"), "ROIs");
  AddLayers(wnd->GetLayers("PointSet"), "Point Sets");
  AddLayers(wnd->GetLayers("CMAT"), "CMAT");
}

void PanelAllLayers::AddLayers(QList<Layer *> layers, const QString &cat_name)
{
  ui->treeWidgetLayers->blockSignals(true);
  if (!layers.isEmpty())
  {
    QTreeWidgetItem* topItem = new QTreeWidgetItem(ui->treeWidgetLayers);
    topItem->setText(0, cat_name);
    for (int i = 0; i < layers.size(); i++)
    {
      QTreeWidgetItem* item = new QTreeWidgetItem(topItem);
      item->setText(0, layers[i]->GetName());
      item->setData(0, Qt::UserRole, QVariant::fromValue((QObject*)layers[i]));
      item->setCheckState(0, layers[i]->IsVisible() ? Qt::Checked : Qt::Unchecked);
    }
    topItem->setExpanded(true);
  }
  ui->treeWidgetLayers->blockSignals(false);
}

void PanelAllLayers::OnLayerRemoved(Layer * removed_ayer)
{
  ui->treeWidgetLayers->blockSignals(true);
  for (int i = 0; i < ui->treeWidgetLayers->topLevelItemCount(); i++)
  {
    QTreeWidgetItem* topItem = ui->treeWidgetLayers->topLevelItem(i);
    for (int j = 0; j < topItem->childCount(); j++)
    {
      QTreeWidgetItem* item = topItem->child(j);
      Layer* layer = qobject_cast<Layer*>(item->data(0, Qt::UserRole).value<QObject*>());
      if (layer == removed_ayer)
      {
        topItem->removeChild(item);
        delete item;
        ui->treeWidgetLayers->blockSignals(false);
        return;
      }
    }
  }
  ui->treeWidgetLayers->blockSignals(false);
}

void PanelAllLayers::OnCurrentItemChanged(QTreeWidgetItem *item)
{
  if (!item)
    return;

  Layer* layer = qobject_cast<Layer*>(item->data(0, Qt::UserRole).value<QObject*>());
  if (!layer)
    return;

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if (layer->IsTypeOf("MRI"))
    mainwnd->GetLayerCollection("MRI")->SetActiveLayer(layer);
  else if (layer->IsTypeOf("Surface"))
    mainwnd->GetLayerCollection("Surface")->SetActiveLayer(layer);
  else if (layer->IsTypeOf("ROI"))
    mainwnd->GetLayerCollection("ROI")->SetActiveLayer(layer);
  else if (layer->IsTypeOf("PointSet"))
    mainwnd->GetLayerCollection("PointSet")->SetActiveLayer(layer);
  else if (layer->IsTypeOf("CMAT"))
    mainwnd->GetLayerCollection("CMAT")->SetActiveLayer(layer);
}

void PanelAllLayers::OnItemChanged(QTreeWidgetItem *item)
{
  Layer* layer = qobject_cast<Layer*>(item->data( 0, Qt::UserRole ).value<QObject*>());
  if ( layer )
  {
  //  layer->SetName( item->text(0) );
    layer->SetVisible( item->checkState( 0 ) == Qt::Checked );
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
  if (!layer)
    return;

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
}
