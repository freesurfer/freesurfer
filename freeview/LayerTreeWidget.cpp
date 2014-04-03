/**
 * @file  LayerTreeWidget.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2014/04/03 20:23:26 $
 *    $Revision: 1.12 $
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
#include <QPainter>
#include <QContextMenuEvent>
#include <QMenu>
#include <QDebug>
#include "MainWindow.h"
#include "ui_MainWindow.h"

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
}

void LayerTreeWidget::drawRow( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const
{
  QTreeWidget::drawRow( painter, option, index );

  Layer* layer = qobject_cast<Layer*>( index.data( Qt::UserRole ).value<QObject*>() );
  if ( layer && layer->IsLocked() )
  {
    QImage img( ":resource/icons/volume_lock.png");
    QRect rc = option.rect;
    rc.setLeft( rc.right() - 20 );
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
      return;
    }
    viewOptions();
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
      return;
    }
    QTreeWidget::mouseReleaseEvent(event);
  }
}

void LayerTreeWidget::contextMenuEvent(QContextMenuEvent *e)
{
  QList<QTreeWidgetItem*> items = selectedItems();
  QList<Layer*> layers;
  foreach (QTreeWidgetItem* item, items)
  {
     Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
     if (layer)
       layers << layer;
  }

  QString type;
  QTreeWidgetItem* item = itemAt(e->pos());
  Layer* layer = NULL;
  if (item)
  {
    layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
    if (layer)
      type = layer->GetPrimaryType();
    else
    {
      type = item->data(0, Qt::UserRole).toString();
    }
  }

  MainWindow* wnd = MainWindow::GetMainWindow();
  QMenu* menu = new QMenu(this);

  if (type == "MRI" || type.isEmpty())
  {
    menu->addAction(wnd->ui->actionNewVolume);
    menu->addAction(wnd->ui->actionLoadVolume);
    if (type == "MRI")
    {
      QAction* act = new QAction("Save All Checked Volumes", this);
      connect(act, SIGNAL(triggered()), this, SLOT(OnSaveVisibleVolumes()));
      menu->addAction(act);
    }
    menu->addSeparator();
  }
  if (type == "Surface" || type.isEmpty())
  {
    menu->addAction(wnd->ui->actionLoadSurface);
    menu->addSeparator();
  }
  if (type == "ROI" || type.isEmpty())
  {
    menu->addAction(wnd->ui->actionNewROI);
    menu->addAction(wnd->ui->actionLoadROI);
    menu->addAction(wnd->ui->actionGoToROI);
    menu->addSeparator();
  }
  if (type == "PointSet" || type.isEmpty())
  {
    menu->addAction(wnd->ui->actionNewPointSet);
    menu->addAction(wnd->ui->actionLoadPointSet);
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
    menu->addSeparator();
    act = new QAction(layers.size() > 1 ? "Show All in Info Panel" : "Show Info", this );
    connect(act, SIGNAL(triggered()), this, SLOT(OnShowAllInfo()));
    menu->addAction(act);
    act = new QAction(layers.size() > 1 ? "Hide All in Info Panel" : "Hide Info", this );
    connect(act, SIGNAL(triggered()), this, SLOT(OnHideAllInfo()));
    menu->addAction(act);

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
    }
  }

  menu->exec(e->globalPos());
}

void LayerTreeWidget::OnShowAll()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
    if (layer)
      layer->Show();
  }
}

void LayerTreeWidget::OnHideAll()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
    if (layer)
      layer->Hide();
  }
}

void LayerTreeWidget::OnLockAll()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
    if (layer)
      layer->Lock(true);
  }
}

void LayerTreeWidget::OnUnlockAll()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
    if (layer)
      layer->Lock(false);
  }
}

void LayerTreeWidget::OnShowAllInfo()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
    if (layer)
      layer->GetProperty()->SetShowInfo(true);
  }
}

void LayerTreeWidget::OnHideAllInfo()
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
    if (layer)
      layer->GetProperty()->SetShowInfo(false);
  }
}

void LayerTreeWidget::OnSetColorMap()
{
  QAction* act = qobject_cast<QAction*>(sender());
  QList<QTreeWidgetItem*> items = this->selectedItems();
  foreach (QTreeWidgetItem* item, items)
  {
    LayerMRI* layer = qobject_cast<LayerMRI*>( item->data(0, Qt::UserRole ).value<QObject*>() );
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
