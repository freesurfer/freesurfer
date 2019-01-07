#ifndef LABELTREEWIDGET_H
#define LABELTREEWIDGET_H

#include <QTreeWidget>

class QTreeWidgetItem;

class LabelTreeWidget : public QTreeWidget
{
  Q_OBJECT
public:
  explicit LabelTreeWidget(QWidget *parent = 0);

  void contextMenuEvent(QContextMenuEvent *e);

  virtual void dropEvent(QDropEvent * event);
  virtual void dragEnterEvent(QDragEnterEvent *event);

signals:
  void MenuGoToCentroid();
  void MenuResample();
  void MenuMoreOps();
  void MenuSaveAs();

public slots:

private:
  QTreeWidgetItem* draggedItem;
};

#endif // LABELTREEWIDGET_H
