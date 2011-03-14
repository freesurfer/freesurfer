#include "LayerTreeWidget.h"
#include "Layer.h"
#include <QPainter>
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
    painter->drawImage( rc.topLeft()+QPoint(0,1), img.scaled( 16, 16) );
  }
}
