#include "SplineTreeWidget.h"
#include <QContextMenuEvent>
#include <QMenu>
#include <QDebug>
#include "SurfaceSpline.h"
#include <QPainter>

SplineTreeWidget::SplineTreeWidget(QWidget *parent) :
  QTreeWidget(parent)
{
}


void SplineTreeWidget::contextMenuEvent(QContextMenuEvent *e)
{
  QTreeWidgetItem* item = currentItem();
  if (item)
  {
    SurfaceSpline* spline = reinterpret_cast<SurfaceSpline*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (spline)
    {
      bool bLocked = spline->IsLocked();
      QMenu* menu = new QMenu(this);
      QAction* act = new QAction(bLocked?"Unlock":"Lock", this);
      connect(act, SIGNAL(triggered()), this, SLOT(OnItemLock()));
      menu->addAction(act);
      menu->exec(e->globalPos());
    }
  }
}

void SplineTreeWidget::OnItemLock()
{
  QTreeWidgetItem* item = currentItem();
  if (item)
  {
    SurfaceSpline* spline = reinterpret_cast<SurfaceSpline*>( item->data(0, Qt::UserRole ).value<quintptr>() );
    if (spline)
      spline->SetLocked(!spline->IsLocked());
    update();
  }
}

void SplineTreeWidget::drawRow( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const
{
  QTreeWidget::drawRow( painter, option, index );

  SurfaceSpline* spline = reinterpret_cast<SurfaceSpline*>( index.data( Qt::UserRole ).value<quintptr>() );
  if ( spline && spline->IsLocked() )
  {
    QImage img( ":resource/icons/volume_lock.png");
    QRect rc = option.rect;
    rc.setLeft( rc.right() - 20 );
    int nsize = qMin(16, rc.height());
    painter->drawImage( rc.topLeft(),
                        img.scaled( nsize, nsize, Qt::KeepAspectRatio, Qt::SmoothTransformation) );
  }
}
