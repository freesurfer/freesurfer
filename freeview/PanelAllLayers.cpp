#include "PanelAllLayers.h"
#include "ui_PanelAllLayers.h"
#include "Layer.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerSurface.h"
#include "LayerPropertySurface.h"
#include "LayerROI.h"
#include "LayerPropertyROI.h"
#include "LayerPointSet.h"
#include "LayerPropertyPointSet.h"
#include <QDebug>
#include <QSettings>

PanelAllLayers::PanelAllLayers(QWidget *parent) :
  QScrollArea(parent),
  ui(new Ui::PanelAllLayers)
{
  ui->setupUi(this);
  MainWindow* wnd = MainWindow::GetMainWindow();
  QStringList layer_types;
  layer_types << "MRI" << "Surface" << "ROI" << "PointSet" << "CMAT" << "Tract" << "FCD";
  foreach (QString type, layer_types)
  {
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerAdded(Layer*)), this, SLOT(OnLayerAdded(Layer*)));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerRemoved(Layer*)), this, SLOT(OnLayerRemoved(Layer*)));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerMoved(Layer*)), this, SLOT(RefreshLayerList()));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerCycled(Layer*)), this, SLOT(RefreshLayerList()));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayersReordered()), this, SLOT(RefreshLayerList()));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerNameChanged()), this, SLOT(OnLayerChanged()));
    connect(wnd->GetLayerCollection(type), SIGNAL(LayerVisibilityChanged()), this, SLOT(OnLayerChanged()));
    connect(wnd->GetLayerCollection(type), SIGNAL(ActiveLayerChanged(Layer*)), this, SLOT(OnActiveLayerChanged(Layer*)));
  }
  ui->treeWidgetLayers->setEditTriggers( QTreeWidget::SelectedClicked );
  ui->stackedWidget->hide();
  connect(ui->treeWidgetLayers, SIGNAL(ToReorderLayers(QList<Layer*>)), this, SIGNAL(ToReorderLayers(QList<Layer*>)));
  connect(wnd, SIGNAL(LinkVolumeRequested(LayerMRI*)), ui->treeWidgetLayers, SLOT(LinkVolume(LayerMRI*)));

  for (int i = 0; i < ui->stackedWidget->count(); i++)
  {
    PanelLayer* panel = qobject_cast<PanelLayer*>(ui->stackedWidget->widget(i));
    if (panel)
      panel->InitializeLayerTreeWidget(ui->treeWidgetLayers);
  }

  QSettings s;
  QVariant v = s.value("ControlPanel/SplitterState");
  if (v.isValid())
  {
    ui->splitterControlPanel->restoreState( v.toByteArray());
  }
  else
  {
    QList<int> sizes;
    sizes << 120 << ui->splitterControlPanel->size().height()-120;
    ui->splitterControlPanel->setSizes(sizes);
  }
}

PanelAllLayers::~PanelAllLayers()
{
  QSettings s;
  s.setValue("ControlPanel/SplitterState", ui->splitterControlPanel->saveState());
  delete ui;
}

void PanelAllLayers::OnActiveLayerChanged(Layer *curlayer)
{
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  if (item && curlayer)
  {
    QObject* sel = reinterpret_cast<QObject*>(item->data(0, Qt::UserRole).value<quintptr>());
    QTreeWidgetItem* topItem = NULL;
    for (int i = 0; i < ui->treeWidgetLayers->topLevelItemCount(); i++)
    {
      item = ui->treeWidgetLayers->topLevelItem(i);
      if (item->data(0, Qt::UserRole).toString() == curlayer->GetPrimaryType())
      {
        topItem = item;
        break;
      }
    }
    if (topItem)
    {
      for (int j = 0; j < topItem->childCount(); j++)
      {
        item = topItem->child(j);
        Layer* layer = reinterpret_cast<Layer*>(item->data(0, Qt::UserRole).value<quintptr>());
        QFont fnt = item->font(0);
        if (sel != curlayer && layer == curlayer)
        {
          ui->treeWidgetLayers->setCurrentItem(item);
        }
        fnt.setBold(layer == curlayer);
        item->setFont(0, fnt);
      }
    }
  }
}

void PanelAllLayers::OnLayerAdded(Layer *added_layer)
{
  QList<Layer*> layers;
  layers << added_layer;
  RefreshLayerList(layers, added_layer);
}

void PanelAllLayers::RefreshLayerList(const QList<Layer *>& selectedLayers_in, Layer* curLayer)
{
  Layer* layer = NULL;
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  if (item)
    layer = reinterpret_cast<Layer*>(item->data( 0, Qt::UserRole ).value<quintptr>());
  if (curLayer)
    layer = curLayer;

  QList<Layer*> selectedLayers = selectedLayers_in;
  if (selectedLayers.isEmpty())
  {
    QList<QTreeWidgetItem*> items = ui->treeWidgetLayers->selectedItems();
    foreach (QTreeWidgetItem* item, items)
    {
      Layer* layer_sel = reinterpret_cast<Layer*>(item->data( 0, Qt::UserRole ).value<quintptr>());
      if (layer_sel)
        selectedLayers << layer_sel;
    }
  }

  ui->treeWidgetLayers->blockSignals(true);
  ui->treeWidgetLayers->clear();
  ui->treeWidgetLayers->blockSignals(false);
  MainWindow* wnd = MainWindow::GetMainWindow();
  AddLayers(wnd->GetLayers("MRI"), "Volumes", wnd->GetActiveLayer("MRI"), selectedLayers, layer );
  AddLayers(wnd->GetLayers("Surface"), "Surfaces", wnd->GetActiveLayer("Surface"), selectedLayers, layer);
  AddLayers(wnd->GetLayers("ROI"), "ROIs", wnd->GetActiveLayer("ROI"), selectedLayers, layer);
  AddLayers(wnd->GetLayers("PointSet"), "Point Sets", wnd->GetActiveLayer("PointSet"), selectedLayers, layer);
  AddLayers(wnd->GetLayers("CMAT"), "CMAT", wnd->GetActiveLayer("CMAT"), selectedLayers, layer);
  AddLayers(wnd->GetLayers("Tract"), "Tract", wnd->GetActiveLayer("Tract"), selectedLayers, layer);
  AddLayers(wnd->GetLayers("FCD"), "FCD", wnd->GetActiveLayer("FCD"), selectedLayers, layer);
}

void PanelAllLayers::AddLayers(QList<Layer *> layers, const QString &cat_name, Layer* activeLayer,
                               const QList<Layer *>& selectedLayers, Layer* curLayer)
{
  ui->treeWidgetLayers->blockSignals(true);
  QTreeWidgetItem* currentItem = NULL;
  QList<QTreeWidgetItem*> selectedItems;
  if (!layers.isEmpty())
  {
    QTreeWidgetItem* topItem = new QTreeWidgetItem(ui->treeWidgetLayers);
    topItem->setText(0, cat_name);
    topItem->setData(0, Qt::UserRole, layers[0]->GetPrimaryType());
    topItem->setFlags(topItem->flags() & (~Qt::ItemIsSelectable));
    for (int i = 0; i < layers.size(); i++)
    {
      QTreeWidgetItem* item = new QTreeWidgetItem();
      item->setText(0, layers[i]->GetName());
      if (layers[i]->IsTypeOf("Surface"))
      {
        LayerSurface* surf = (LayerSurface*)layers[i];
        SetItemColor(item, surf->GetProperty()->GetEdgeColor());
        connect(surf->GetProperty(), SIGNAL(EdgeColorChanged()), this, SLOT(OnLayerChanged()), Qt::UniqueConnection);
      }
      else if (layers[i]->IsTypeOf("ROI"))
      {
        LayerROI* roi = (LayerROI*)layers[i];
        SetItemColor(item, roi->GetProperty()->GetColor());
        connect(roi->GetProperty(), SIGNAL(ColorMapChanged()), this, SLOT(OnLayerChanged()), Qt::UniqueConnection);
      }
      else if (layers[i]->IsTypeOf("PointSet"))
      {
        LayerPointSet* layer = (LayerPointSet*)layers[i];
        SetItemColor(item, layer->GetProperty()->GetColor());
        connect(layer->GetProperty(), SIGNAL(ColorChanged()), this, SLOT(OnLayerChanged()), Qt::UniqueConnection);
      }
      item->setData(0, Qt::UserRole, QVariant::fromValue(reinterpret_cast<quintptr>(layers[i])));
      item->setCheckState(0, layers[i]->IsVisible() ? Qt::Checked : Qt::Unchecked);
      item->setFlags( item->flags() | Qt::ItemIsEditable | Qt::ItemIsDragEnabled);
      topItem->addChild(item);
      if (layers[i] == curLayer)
      {
        currentItem = item;
      }
      if (layers[i] == activeLayer)
      {
        QFont fnt = item->font(0);
        fnt.setBold(true);
        item->setFont(0, fnt);
      }
      if (selectedLayers.contains(layers[i]))
      {
        selectedItems << item;
      }
      if (layers[i]->GetAboutToDelete())
        item->setHidden(true);
    }
    topItem->setExpanded(true);
  }
  ui->treeWidgetLayers->blockSignals(false);
  if (currentItem)
    ui->treeWidgetLayers->setCurrentItem(currentItem);
  if (!selectedItems.isEmpty())
  {
    foreach (QTreeWidgetItem* item, selectedItems)
      item->setSelected(true);
  }
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
      if ( reinterpret_cast<Layer*>(item->data(0, Qt::UserRole).value<quintptr>()) == removed_layer)
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

  Layer* layer = reinterpret_cast<Layer*>(item->data(0, Qt::UserRole).value<quintptr>());
  QString type;
  if (layer)
  {
    MainWindow* mainwnd = MainWindow::GetMainWindow();
    type = layer->GetPrimaryType();
    LayerCollection* lc = mainwnd->GetLayerCollection(type);
    if (lc)
    {
      lc->SetActiveLayer(layer);
    }
  }
  else
  {
    type = item->data(0, Qt::UserRole).toString();
  }

  PanelLayer* panel = SetCurrentPanel(type);
  if (panel)
  {
    panel->SetCurrentLayer(layer);
    emit CurrentLayerSelected(layer);
  }
}

void PanelAllLayers::OnItemChanged(QTreeWidgetItem *item)
{
  Layer* layer = reinterpret_cast<Layer*>(item->data( 0, Qt::UserRole ).value<quintptr>());
  if ( layer )
  {
    if (item->text(0) != layer->GetName())
      layer->SetName( item->text(0) );
    if (layer->IsVisible() != (item->checkState( 0 ) == Qt::Checked ))
      layer->SetVisible( item->checkState( 0 ) == Qt::Checked );
  }
}

void PanelAllLayers::OnItemSelectionChanged()
{
  // first disconnect all previous layer connections
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  if (item)
  {
    Layer* layer = reinterpret_cast<Layer*>(item->data(0, Qt::UserRole).value<quintptr>());
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
      {
        panel->ConnectLayer(l);
      }
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
      Layer* layer = reinterpret_cast<Layer*>(item->data(0, Qt::UserRole).value<quintptr>());
      if (layer)
      {
        item->setText(0, layer->GetName());
        item->setCheckState(0, layer->IsVisible() ?  Qt::Checked : Qt::Unchecked);
        if (layer->IsTypeOf("Surface"))
        {
          LayerSurface* surf = (LayerSurface*)layer;
          SetItemColor(item, surf->GetProperty()->GetEdgeColor());
          connect(surf->GetProperty(), SIGNAL(EdgeColorChanged()), this, SLOT(OnLayerChanged()), Qt::UniqueConnection);
        }
        else if (layer->IsTypeOf("ROI"))
        {
          LayerROI* roi = (LayerROI*)layer;
          SetItemColor(item, roi->GetProperty()->GetColor());
          connect(roi->GetProperty(), SIGNAL(ColorMapChanged()), this, SLOT(OnLayerChanged()), Qt::UniqueConnection);
        }
        else if (layer->IsTypeOf("PointSet"))
        {
          LayerPointSet* ps = (LayerPointSet*)layer;
          SetItemColor(item, ps->GetProperty()->GetColor());
          connect(ps->GetProperty(), SIGNAL(ColorChanged()), this, SLOT(OnLayerChanged()), Qt::UniqueConnection);
        }
      }
    }
  }
  ui->treeWidgetLayers->blockSignals(false);
}

void PanelAllLayers::SetItemColor(QTreeWidgetItem *item, double *rgb)
{
  QPixmap pix(13, 13);
  pix.fill( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
  item->setIcon(0, QIcon(pix) );
}

void PanelAllLayers::OnItemDoubleClicked(QTreeWidgetItem *item)
{
  Layer* layer = reinterpret_cast<Layer*>(item->data( 0, Qt::UserRole ).value<quintptr>());
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
    Layer* layer = reinterpret_cast<Layer*>(item->data( 0, Qt::UserRole ).value<quintptr>());
    if (layer)
    {
      QStringList layer_types;
      layer_types << "MRI" << "Surface" << "ROI" << "PointSet" << "CMAT" << "Tract" << "FCD";
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
    Layer* layer = reinterpret_cast<Layer*>(item->data( 0, Qt::UserRole ).value<quintptr>());
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

void PanelAllLayers::SelectAllLayers()
{
  ui->treeWidgetLayers->SelectAll();
}

void PanelAllLayers::DeselectAllLayers()
{
  ui->treeWidgetLayers->DeselectAll();
}

void PanelAllLayers::SetSelectedLayers(const QList<int> &layer_ids)
{
  ui->treeWidgetLayers->SetSelectedLayers(layer_ids);
}

PanelLayer* PanelAllLayers::GetPanel(const QString &layer_type)
{
  for (int i = 0; i < ui->stackedWidget->count(); i++)
  {
    PanelLayer* panel = qobject_cast<PanelLayer*>(ui->stackedWidget->widget(i));
    if (panel && panel->GetLayerType() == layer_type)
    {
      ui->stackedWidget->setCurrentWidget(panel);
      return panel;
    }
  }
  return NULL;
}
