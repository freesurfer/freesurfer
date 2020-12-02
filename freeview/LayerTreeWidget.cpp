/*
 * Original Author: Ruopeng Wang
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
#include "LayerTreeWidget.h"
#include "Layer.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerSurface.h"
#include "LayerPropertySurface.h"
#include <QPainter>
#include <QContextMenuEvent>
#include <QMenu>
#include <QDebug>
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QDropEvent>

QRect MyItemDelegate::GetCheckBoxRect(const QModelIndex &index, const QStyleOptionViewItem& option) const
{
  QRect CheckBox = rect (option, index, Qt::CheckStateRole);
  QRect Icon = rect (option, index, Qt::DecorationRole) ;
  QRect Text = rect (option, index, Qt::DisplayRole) ;

  doLayout (option, &CheckBox, &Icon, &Text, true) ;

  QRect VisualRect = ParentView -> visualRect (index) ;
  CheckBox.translate (VisualRect.topLeft()) ;
  Icon.translate (VisualRect.topLeft()) ;
  Text.translate (VisualRect.topLeft()) ;
  return CheckBox;
}

LayerTreeWidget::LayerTreeWidget(QWidget *parent) :
  QTreeWidget(parent)
{
  m_itemDelegate = new MyItemDelegate(this);
  setItemDelegate(m_itemDelegate);

  //  QAction* act = new QAction("Select All", this);
  //  act->setShortcut(QKeySequence("Ctrl+A"));
  //  act->setShortcutContext(Qt::WidgetWithChildrenShortcut);
  //  connect(act, SIGNAL(triggered()), SLOT(selectAll()));
  //  this->addAction(act);

  setMouseTracking(true);
  setDragEnabled(true);
  viewport()->setAcceptDrops(true);
  setDropIndicatorShown(true);
  setDragDropMode(QAbstractItemView::InternalMove);
}

void LayerTreeWidget::drawRow( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const
{
  QTreeWidget::drawRow( painter, option, index );

  Layer* layer = reinterpret_cast<Layer*>( index.data( Qt::UserRole ).value<quintptr>() );
  QRect rc = option.rect;
  rc.setLeft( rc.right() - 20 );
  QTreeWidgetItem* item = itemAt(rc.center());
  if (item)
      item->setData(0, Qt::UserRole+10, rc);

  if ( layer && layer->IsLocked())
  {
    QImage img(":resource/icons/volume_lock.png");
    int nsize = qMin(16, rc.height());
    painter->drawImage( rc.topLeft(),
                        img.scaled( nsize, nsize, Qt::KeepAspectRatio, Qt::SmoothTransformation) );
  }
}

void LayerTreeWidget::ForceUpdate()
{
  this->setDirtyRegion(QRegion(rect()));
  update();
}

void LayerTreeWidget::mousePressEvent(QMouseEvent *event)
{
  if (event->button()== Qt::RightButton)
    return;
  else
  {
    QTreeWidgetItem* item = itemAt(event->pos());
    if (rectCheckbox.isEmpty() || rectCheckbox.width() <= 0)
    {
      QStyleOptionViewItem option = viewOptions();
      if (item)
        rectCheckbox = m_itemDelegate->GetCheckBoxRect(indexFromItem(item), option);
    }

    if (item && item->childCount() == 0
        && event->x() < rectCheckbox.right() && event->x() > rectCheckbox.left())
    {
      m_bCheckBoxClicked = true;
      return;
    }

    bool bClickToLock = MainWindow::GetMainWindow()->GetSetting("ClickToLock").toBool();
    Layer* layer = NULL;
    if (item)
      layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );

    if ( layer && (layer->IsLocked() || bClickToLock) && item->data(0, Qt::UserRole+10).toRect().contains(event->pos()))
      return;

    QTreeWidget::mousePressEvent(event);
  }
}

void LayerTreeWidget::mouseReleaseEvent(QMouseEvent *event)
{
  if (event->button() == Qt::LeftButton)
  {
    QTreeWidgetItem* item = itemAt(event->pos());
    if (item && item->childCount() == 0
        && event->x() < rectCheckbox.right() && event->x() > rectCheckbox.left())
    {
      item->setCheckState(0, item->checkState(0) == Qt::Checked ? Qt::Unchecked : Qt::Checked);
      m_bCheckBoxClicked = false;
      return;
    }

    bool bClickToLock = MainWindow::GetMainWindow()->GetSetting("ClickToLock").toBool();
    Layer* layer = NULL;
    if (item)
      layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );

    if ( layer && (layer->IsLocked() || bClickToLock) && item->data(0, Qt::UserRole+10).toRect().contains(event->pos()))
    {
      layer->Lock(!layer->IsLocked());
      return;
    }

    QTreeWidget::mouseReleaseEvent(event);
  }
  m_bCheckBoxClicked = false;
}

void LayerTreeWidget::mouseMoveEvent(QMouseEvent *event)
{
  if (m_bCheckBoxClicked)
    return;

  QTreeWidgetItem* item = itemAt(event->pos());
  Layer* layer = NULL;
  if (item)
    layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );

  bool bClickToLock = MainWindow::GetMainWindow()->GetSetting("ClickToLock").toBool();
  if ( layer && (layer->IsLocked() || bClickToLock) && item->data(0, Qt::UserRole+10).toRect().contains(event->pos()) )
  {
    setCursor(Qt::PointingHandCursor);
    return;
  }
  else
    unsetCursor();

  QTreeWidget::mouseMoveEvent(event);
}

void LayerTreeWidget::contextMenuEvent(QContextMenuEvent *e)
{
  QList<QTreeWidgetItem*> items = selectedItems();
  QList<Layer*> layers;
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
      layers << layer;
  }

  QString type;
  QTreeWidgetItem* item = itemAt(e->pos());
  Layer* layer = NULL;
  if (item)
  {
    layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
      type = layer->GetPrimaryType();
    else
    {
      type = item->data(0, Qt::UserRole).toString();
    }
  }

  MainWindow* wnd = MainWindow::GetMainWindow();
  QMenu* menu = new QMenu(this);

  //  if (layer)
  //  {
  //    QAction* act = new QAction("Rename", this);
  //    connect(act, SIGNAL(triggered()), this, SLOT(OnEditName()));
  //    menu->addAction(act);
  //    menu->addSeparator();
  //  }

  if (type == "MRI" || type == "Surface")
  {
    if (type == "MRI")
      wnd->ui->actionViewLayerInfo->setText("View Volume Info...");
    else if (type == "Surface")
      wnd->ui->actionViewLayerInfo->setText("View Surface Info...");
    menu->addAction(wnd->ui->actionViewLayerInfo);
    menu->addSeparator();
  }
  if (type == "MRI" || type.isEmpty())
  {
    menu->addAction(wnd->ui->actionNewVolume);
    menu->addAction(wnd->ui->actionLoadVolume);
    menu->addAction(wnd->ui->actionReloadVolume);
    if (type == "MRI")
    {
      QAction* act = new QAction("Save All Checked Volumes", this);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSaveVisibleVolumes()));
      menu->addAction(act);
      menu->addSeparator();
      menu->addAction(wnd->ui->actionLoadTransform);
      if (layer && ((LayerMRI*)layer)->HasReg())
        menu->addAction(wnd->ui->actionUnloadTransform);
    }
    menu->addSeparator();
  }
  if (type == "Surface" || type.isEmpty())
  {
    menu->addAction(wnd->ui->actionLoadSurface);
    menu->addAction(wnd->ui->actionReloadSurface);
    menu->addSeparator();
    if (type == "Surface")
    {
      menu->addAction(wnd->ui->actionLoadPatch);
      menu->addAction(wnd->ui->actionSavePatchAs);
      menu->addSeparator();
    }
  }
  if (type == "ROI" || type.isEmpty())
  {
    menu->addAction(wnd->ui->actionNewROI);
    menu->addAction(wnd->ui->actionLoadROI);
    menu->addAction(wnd->ui->actionReloadROI);
    menu->addSeparator();
    menu->addAction(wnd->ui->actionGoToROI);
    menu->addAction(wnd->ui->actionShowLabelStats);
    menu->addSeparator();
  }
  if (type == "PointSet" || type.isEmpty())
  {
    menu->addAction(wnd->ui->actionNewPointSet);
    menu->addAction(wnd->ui->actionLoadPointSet);
    menu->addAction(wnd->ui->actionReloadPointSet);
    menu->addSeparator();
  }

  if (layers.size() > 0 && items.contains(item))
  {
    QAction* act = new QAction(layers.size() > 1 ? "Show All" : "Show", this );
    connect(act, SIGNAL(triggered()), this, SLOT(OnShowAll()));
    menu->addAction(act);
    act = new QAction(layers.size() > 1 ? "Hide All" : "Hide", this );
    connect(act, SIGNAL(triggered()), this, SLOT(OnHideAll()));
    menu->addAction(act);
    menu->addSeparator();
    act = new QAction(layers.size() > 1 ? "Lock All" : "Lock", this );
    connect(act, SIGNAL(triggered()), this, SLOT(OnLockAll()));
    menu->addAction(act);
    act = new QAction(layers.size() > 1 ? "Unlock All" : "Unlock", this );
    connect(act, SIGNAL(triggered()), this, SLOT(OnUnlockAll()));
    menu->addAction(act);

    act = new QAction("Lock Others", this );
    connect(act, SIGNAL(triggered()), this, SLOT(OnLockOthers()));
    menu->addAction(act);
    act = new QAction("Unlock Others", this );
    connect(act, SIGNAL(triggered()), this, SLOT(OnUnlockOthers()));
    menu->addAction(act);

    if (layers[0]->IsTypeOf("MRI"))
    {
        menu->addSeparator();
        if (layers.size() > 1)
        {
            act = new QAction("Link Volumes", this );
            connect(act, SIGNAL(triggered()), this, SLOT(OnLinkVolumes()));
            menu->addAction(act);
        }
        act = new QAction("Unlink Volumes", this );
        connect(act, SIGNAL(triggered()), this, SLOT(OnUnlinkVolumes()));
        menu->addAction(act);
    }

    if (layers[0]->IsTypeOf("MRI") || layers[0]->IsTypeOf("Surface"))
    {
      menu->addSeparator();
      act = new QAction(layers.size() > 1 ? "Show All in Info Panel" : "Show Info", this );
      connect(act, SIGNAL(triggered()), this, SLOT(OnShowAllInfo()));
      menu->addAction(act);
      act = new QAction(layers.size() > 1 ? "Hide All in Info Panel" : "Hide Info", this );
      connect(act, SIGNAL(triggered()), this, SLOT(OnHideAllInfo()));
      menu->addAction(act);
    }

    if (layers[0]->GetEndType() == "MRI")
    {
      menu->addSeparator();
      QMenu* submenu = new QMenu("Color Map", this);
      menu->addMenu(submenu);
      int nColorMap = ((LayerMRI*)layers[0])->GetProperty()->GetColorMap();
      act = new QAction("Grayscale", this);
      act->setData(LayerPropertyMRI::Grayscale);
      act->setCheckable(true);
      act->setChecked(nColorMap == LayerPropertyMRI::Grayscale);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
      submenu->addAction(act);
      act = new QAction("Lookup Table", this);
      act->setData(LayerPropertyMRI::LUT);
      act->setCheckable(true);
      act->setChecked(nColorMap == LayerPropertyMRI::LUT);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
      submenu->addAction(act);
      act = new QAction("Heat", this);
      act->setData(LayerPropertyMRI::Heat);
      act->setCheckable(true);
      act->setChecked(nColorMap == LayerPropertyMRI::Heat);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
      submenu->addAction(act);
      act = new QAction("Jet", this);
      act->setData(LayerPropertyMRI::Jet);
      act->setCheckable(true);
      act->setChecked(nColorMap == LayerPropertyMRI::Jet);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
      submenu->addAction(act);
      act = new QAction("GE Color", this);
      act->setData(LayerPropertyMRI::GEColor);
      act->setCheckable(true);
      act->setChecked(nColorMap == LayerPropertyMRI::GEColor);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
      submenu->addAction(act);
      act = new QAction("NIH", this);
      act->setData(LayerPropertyMRI::NIH);
      act->setCheckable(true);
      act->setChecked(nColorMap == LayerPropertyMRI::NIH);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
      submenu->addAction(act);
      act = new QAction("PET", this);
      act->setData(LayerPropertyMRI::PET);
      act->setCheckable(true);
      act->setChecked(nColorMap == LayerPropertyMRI::PET);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
      submenu->addAction(act);

      if (((LayerMRI*)layers[0])->GetRefVolume())
      {
        act = new QAction("Apply Transformation...", this);
        connect(act, SIGNAL(triggered(bool)), MainWindow::GetMainWindow(), SLOT(OnApplyVolumeTransform()));
        menu->addAction(act);
      }

      if (((LayerMRI*)layers[0])->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT)
      {
          act = new QAction("Export label stats...", this);
          connect(act, SIGNAL(triggered(bool)), MainWindow::GetMainWindow(), SLOT(OnExportLabelStats()));
          menu->addSeparator();
          menu->addAction(act);
      }
    }
    else if (layers[0]->GetEndType() == "Surface")
    {
      if (((LayerSurface*)layers[0])->IsContralateralPossible())
      {
        menu->addSeparator();
        act = new QAction("Go To Contralateral Point", this);
        connect(act, SIGNAL(triggered(bool)), MainWindow::GetMainWindow(), SLOT(GoToContralateralPoint()));
        menu->addAction(act);
      }
    }
    else if (layers[0]->GetEndType() == "PointSet")
    {
      menu->addSeparator();
      act = new QAction("Go to Centroid", this);
      connect(act, SIGNAL(triggered()), MainWindow::GetMainWindow(), SLOT(OnGoToPointSet()));
      menu->addAction(act);
    }
  }

  menu->exec(e->globalPos());
}

void LayerTreeWidget::OnShowAll()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
      layer->Show();
  }
}

void LayerTreeWidget::OnHideAll()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
      layer->Hide();
  }
}

void LayerTreeWidget::OnLockAll()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
      layer->Lock(true);
  }
}

void LayerTreeWidget::OnUnlockAll()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
      layer->Lock(false);
  }
}

void LayerTreeWidget::OnLockOthers()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  QList<Layer*> selected_layers;
  QString type;
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
    {
      layer->Lock(false);
      selected_layers << layer;
      type = layer->GetPrimaryType();
    }
  }
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers(type);
  foreach (Layer* layer, layers)
  {
    if (!selected_layers.contains(layer))
      layer->Lock(true);
  }
}

void LayerTreeWidget::OnUnlockOthers()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  QList<Layer*> selected_layers;
  QString type;
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
    {
      layer->Lock(true);
      selected_layers << layer;
      type = layer->GetPrimaryType();
    }
  }
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers(type);
  foreach (Layer* layer, layers)
  {
    if (!selected_layers.contains(layer))
      layer->Lock(false);
  }
}

void LayerTreeWidget::OnShowAllInfo()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
      layer->GetProperty()->SetShowInfo(true);
  }
}

void LayerTreeWidget::OnHideAllInfo()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
      layer->GetProperty()->SetShowInfo(false);
  }
}

void LayerTreeWidget::OnEditName()
{
  QTreeWidgetItem* item = this->currentItem();
  if (item)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer)
    {

    }
  }
}

void LayerTreeWidget::OnSetColorMap()
{
  QAction* act = qobject_cast<QAction*>(sender());
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    LayerMRI* layer = reinterpret_cast<LayerMRI*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (layer && layer->GetEndType() == "MRI")
      layer->GetProperty()->SetColorMap(act->data().toInt());
  }
}

void LayerTreeWidget::OnSaveVisibleVolumes()
{
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
  QList<Layer*> visibles;
  foreach (Layer* layer, layers)
  {
    if (layer->IsVisible() && ((LayerMRI*)layer)->IsModified())
    {
      visibles << layer;
    }
  }
  MainWindow::GetMainWindow()->SaveLayers(visibles);
}

void LayerTreeWidget::SelectAll()
{
  QList<QTreeWidgetItem*> items;
  for (int i = 0; i < topLevelItemCount(); i++)
  {
    QTreeWidgetItem* topItem = topLevelItem(i);
    for (int j = 0; j < topItem->childCount(); j++)
      items << topItem->child(j);
  }
  foreach (QTreeWidgetItem* item, items)
  {
    item->setSelected(true);
  }
}

void LayerTreeWidget::DeselectAll()
{
  QList<QTreeWidgetItem*> items;
  for (int i = 0; i < topLevelItemCount(); i++)
  {
    QTreeWidgetItem* topItem = topLevelItem(i);
    for (int j = 0; j < topItem->childCount(); j++)
      items << topItem->child(j);
  }
  foreach (QTreeWidgetItem* item, items)
  {
    item->setSelected(false);
  }
}

void LayerTreeWidget::SetSelectedLayers(const QList<int>& layer_ids)
{
  QList<QTreeWidgetItem*> items;
  for (int i = 0; i < topLevelItemCount(); i++)
  {
    QTreeWidgetItem* topItem = topLevelItem(i);
    for (int j = 0; j < topItem->childCount(); j++)
      items << topItem->child(j);
  }
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    item->setSelected(layer && layer_ids.contains(layer->GetID()));
  }
}

bool LayerTreeWidget::event(QEvent *e)
{
  if (e->type() == QEvent::ShortcutOverride || e->type() == QEvent::KeyPress)
  {
    QKeyEvent* ke = static_cast<QKeyEvent*>(e);
    if ( ke )
    {
      if (ke->key()== Qt::Key_A && ke->modifiers() & Qt::ControlModifier)
      {
        SelectAll();
        return true;
      }
    }
  }
  return QTreeWidget::event(e);
}

void LayerTreeWidget::dropEvent(QDropEvent *event)
{
  QModelIndex droppedIndex = indexAt( event->pos() );

  if ( droppedIndex.isValid() )
  {
    DropIndicatorPosition drop_pos = dropIndicatorPosition();
    QTreeWidgetItem* itemTo = itemAt(event->pos());
    QList<QTreeWidgetItem*> items = this->selectedItems(), itemsFrom;
    QTreeWidgetItem* itemCur = this->currentItem();
    QString type;
    if (itemCur && itemCur->parent())
    {
      QTreeWidgetItem* parent = itemCur->parent();
      type = parent->data(0, Qt::UserRole).toString();
      foreach (QTreeWidgetItem* item, items)
      {
        Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
        if (layer && layer->GetPrimaryType() == type)
        {
          itemsFrom << item;
        }
      }
    }
    bool bReordered = false;
    if (!itemsFrom.isEmpty() && itemTo && !itemsFrom.contains(itemTo))
    {
      QTreeWidgetItem* parent = itemTo->parent();
      if (parent)
      {
        if (parent->data(0, Qt::UserRole).toString() == type)
        {
          foreach (QTreeWidgetItem* item, itemsFrom)
          {
            parent->removeChild(item);
          }
          parent->insertChildren(parent->indexOfChild(itemTo) + (drop_pos == QTreeWidget::BelowItem ? 1 : 0), itemsFrom);
        }
      }
      else if (itemTo->data(0, Qt::UserRole).toString() == type)
      {
        foreach (QTreeWidgetItem* item, itemsFrom)
        {
          itemTo->removeChild(item);
        }
        itemTo->insertChildren(0, itemsFrom);
      }
      bReordered = true;
    }
    if (itemCur)
      setCurrentItem(itemCur);

    foreach (QTreeWidgetItem* item, itemsFrom)
    {
      item->setSelected(true);
    }

    if (bReordered)
    {
      if (itemTo->parent())
        itemTo = itemTo->parent();
      QList<Layer*> layers;
      for (int i = 0; i < itemTo->childCount(); i++)
      {
        QTreeWidgetItem* item = itemTo->child(i);
        Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
        if (layer)
        {
          layers << layer;
        }
      }
      emit ToReorderLayers(layers);
    }

    //        QTreeWidgetItem* itemFrom = NULL;
    //        QByteArray encoded = event->mimeData()->data("application/x-qabstractitemmodeldatalist");
    //        QDataStream stream(&encoded, QIODevice::ReadOnly);
    //        if (!stream.atEnd())
    //        {
    //            int row, col;
    //            QMap<int,  QVariant> roleDataMap;
    //            stream >> row >> col >> roleDataMap;
    //            itemFrom = itemFromIndex(rootIndex().child(row, col));
    //        }
    //        qDebug() << itemTo << itemFrom;
  }

  //    QTreeWidget::dropEvent(event);
}


void LayerTreeWidget::OnLinkVolumes()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  m_linkedVolumes.clear();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = reinterpret_cast<Layer*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    LayerMRI* mri = qobject_cast<LayerMRI*>(layer);
    if (mri)
      m_linkedVolumes << mri;
  }
}

void LayerTreeWidget::LinkVolume(LayerMRI *vol)
{
  if (!m_linkedVolumes.contains(vol))
    m_linkedVolumes << vol;
}

void LayerTreeWidget::OnUnlinkVolumes()
{
    m_linkedVolumes.clear();
}
