/**
 * @file  LayerTreeWidget.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
 *    $Revision: 1.4 $
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
