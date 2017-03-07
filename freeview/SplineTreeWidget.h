#ifndef SplineTreeWidget_H
#define SplineTreeWidget_H

#include <QTreeWidget>

class SplineTreeWidget : public QTreeWidget
{
  Q_OBJECT
public:
  explicit SplineTreeWidget(QWidget *parent = 0);

  void contextMenuEvent(QContextMenuEvent *e);

signals:

public slots:
  void OnItemLock();

protected:
  void drawRow ( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const;

};

#endif // SplineTreeWidget_H
