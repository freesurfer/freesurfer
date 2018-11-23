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
    act->setData("go_to_centroid");
    connect(act, SIGNAL(triggered()), this, SLOT(OnMenuTriggered()));
    menu->addAction(act);
    act = new QAction("Resample", this);
    act->setData("resample");
    connect(act, SIGNAL(triggered()), this, SLOT(OnMenuTriggered()));
    menu->addAction(act);
    menu->exec(e->globalPos());
  }
}

void LabelTreeWidget::OnMenuTriggered()
{
  QAction* act = qobject_cast<QAction*>(sender());
  if (act)
  {
    QString act_str = act->data().toString();
    if (act_str == "go_to_centroid")
      emit MenuGoToCentroid();
    else if (act_str == "resample")
      emit MenuResample();
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
