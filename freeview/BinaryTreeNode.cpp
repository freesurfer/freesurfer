#include "BinaryTreeNode.h"
#include "BinaryTreeEdge.h"
#include "BinaryTreeView.h"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsSceneHoverEvent>
#include <QPainter>
#include <QStyleOption>
#include <iostream>

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsPixmapItem>
#include <QDebug>

BinaryTreeNode::BinaryTreeNode(BinaryTreeView *graphWidget)
  : graph(graphWidget)
{
  setFlag(ItemIsSelectable);
  //setFlag(ItemIsMovable);
  setFlag(ItemSendsGeometryChanges);
  setCacheMode(DeviceCoordinateCache);
//  setZValue(-1);
  setAcceptHoverEvents(true);
  setCursor(Qt::PointingHandCursor);
}

void BinaryTreeNode::addEdge(BinaryTreeEdge *edge)
{
  edgeList << edge;
  edge->adjust();
}


QList<BinaryTreeEdge *> BinaryTreeNode::edges() const
{
  return edgeList;
}

bool BinaryTreeNode::advance()
{
  if (newPos == pos())
    return false;

  setPos(newPos);
  return true;
}

QRectF BinaryTreeNode::boundingRect() const
{
  qreal adjust = 2;
  return QRectF( -10 - adjust, -10 - adjust, 23 + adjust, 23 + adjust);
}

QPainterPath BinaryTreeNode::shape() const
{
  QPainterPath path;
  path.addEllipse(-10, -10, 20, 20);
  return path;
}

void BinaryTreeNode::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *)
{
//  painter->setPen(Qt::NoPen);
//  painter->setBrush(Qt::darkGray);
//  painter->drawEllipse(-3.5, -3.5, 5, 5);

  QRadialGradient gradient(-3, -3, 10);
  if (option->state & QStyle::State_Sunken) {
    gradient.setCenter(3, 3);
    gradient.setFocalPoint(3, 3);
    gradient.setColorAt(1, QColor(Qt::yellow).light(120));
    gradient.setColorAt(0, QColor(Qt::darkYellow).light(120));
  } else {
    gradient.setColorAt(0, Qt::yellow);
    gradient.setColorAt(1, Qt::darkYellow);
  }

  painter->setBrush(gradient);

  painter->setPen(QPen(Qt::black, 0));
  //painter->drawEllipse(-10, -10, 20, 20);
  painter->drawEllipse(-5, -5, 10, 10);

}

void BinaryTreeNode::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  QGraphicsItem::mousePressEvent(event);
  update();
}

void BinaryTreeNode::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
  update();
  QGraphicsItem::mouseReleaseEvent(event);
}

