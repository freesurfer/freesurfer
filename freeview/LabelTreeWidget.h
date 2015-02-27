#ifndef LABELTREEWIDGET_H
#define LABELTREEWIDGET_H

#include <QTreeWidget>

class LabelTreeWidget : public QTreeWidget
{
    Q_OBJECT
public:
    explicit LabelTreeWidget(QWidget *parent = 0);

  void contextMenuEvent(QContextMenuEvent *e);

signals:
  void MenuGoToCentroid();

public slots:
  void OnMenuTriggered();

};

#endif // LABELTREEWIDGET_H
