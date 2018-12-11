#include "LabelTreeWidget.h"
#include <QContextMenuEvent>
#include <QMenu>
#include <QDebug>

LabelTreeWidget::LabelTreeWidget(QWidget *parent) :
  QTreeWidget(parent), draggedItem(NULL)
{

}

void LabelTreeWidget::contextMenuEvent(QContextMenuEvent *e)
{
  QTreeWidgetItem* item = currentItem();
  if (item)
  {
    QMenu* menu = new QMenu(this);
    QAction* act = new QAction("Go To Centroid", this);
    connect(act, SIGNAL(triggered()), this, SIGNAL(MenuGoToCentroid()));
    menu->addAction(act);
    act = new QAction("Resample", this);
    connect(act, SIGNAL(triggered()), this, SIGNAL(MenuResample()));
    menu->addAction(act);
    act = new QAction("Dilate/Erode/Open/Close...", this);
    connect(act, SIGNAL(triggered()), this, SIGNAL(MenuMoreOps()));
    menu->addAction(act);
    menu->addSeparator();
    act = new QAction("Save As...", this);
    connect(act, SIGNAL(triggered()), this, SIGNAL(MenuSaveAs()));
    menu->addAction(act);
    menu->exec(e->globalPos());
  }
}

void LabelTreeWidget::dragEnterEvent(QDragEnterEvent *event)
{
  draggedItem = currentItem();
  QTreeWidget::dragEnterEvent(event);
}

void LabelTreeWidget::dropEvent(QDropEvent *event)
{
  QModelIndex droppedIndex = indexAt(event->pos());
  if ( !droppedIndex.isValid() )
    return;

  qDebug() << draggedItem << droppedIndex;
  if (draggedItem)
  {
    QTreeWidgetItem* dParent = draggedItem->parent();
    if (dParent)
    {
      if (itemFromIndex(droppedIndex.parent()) != dParent)
        return;
      dParent->removeChild(draggedItem);
      dParent->insertChild(droppedIndex.row(), draggedItem);
    }
  }
}
