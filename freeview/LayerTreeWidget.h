/**
 * @file  LayerTreeWidget.h
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
#ifndef LAYERTREEWIDGET_H
#define LAYERTREEWIDGET_H

#include <QTreeWidget>

class LayerTreeWidget : public QTreeWidget
{
  Q_OBJECT
public:
  explicit LayerTreeWidget(QWidget *parent = 0);

signals:

public slots:

protected:
  void drawRow ( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const;
};

#endif // LAYERTREEWIDGET_H
