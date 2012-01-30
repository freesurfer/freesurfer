/**
 * @file  LayerTreeWidget.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/01/30 20:57:05 $
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
#include "LayerTreeWidget.h"
#include "Layer.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include <QPainter>
#include <QContextMenuEvent>
#include <QMenu>
#include <QDebug>

LayerTreeWidget::LayerTreeWidget(QWidget *parent) :
  QTreeWidget(parent)
{
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
    QTreeWidget::mousePressEvent(event);
}

void LayerTreeWidget::contextMenuEvent(QContextMenuEvent *e)
{
  QList<QTreeWidgetItem*> items = this->selectedItems();
  QList<Layer*> layers;
  foreach (QTreeWidgetItem* item, items)
  {
    layers << qobject_cast<Layer*>( item->data(0, Qt::UserRole ).value<QObject*>() );
  }
  if (layers.isEmpty())
    return;

  QMenu* menu = new QMenu(this);
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
    act = new QAction("Grayscale", this);
    act->setData(LayerPropertyMRI::Grayscale);
    connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
    submenu->addAction(act);
    act = new QAction("Lookup Table", this);
    act->setData(LayerPropertyMRI::LUT);
    connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
    submenu->addAction(act);
    act = new QAction("Heat", this);
    act->setData(LayerPropertyMRI::Heat);
    connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
    submenu->addAction(act);
    act = new QAction("Jet", this);
    act->setData(LayerPropertyMRI::Jet);
    connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
    submenu->addAction(act);
    act = new QAction("GE Color", this);
    act->setData(LayerPropertyMRI::GEColor);
    connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
    submenu->addAction(act);
    act = new QAction("NIH", this);
    act->setData(LayerPropertyMRI::NIH);
    connect(act, SIGNAL(triggered()), this, SLOT(OnSetColorMap()));
    submenu->addAction(act);
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
