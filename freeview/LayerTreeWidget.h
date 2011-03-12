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
