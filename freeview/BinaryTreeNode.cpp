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

#define NODE_RADIUS 5

BinaryTreeNode::BinaryTreeNode(BinaryTreeView *graphWidget)
  : m_bHighlighted(false)
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
  return QRectF( -NODE_RADIUS - adjust, -NODE_RADIUS - adjust, NODE_RADIUS*2 + adjust, NODE_RADIUS*2 + adjust);
}

QPainterPath BinaryTreeNode::shape() const
{
  QPainterPath path;
  path.addEllipse(-NODE_RADIUS-1, -NODE_RADIUS-1, NODE_RADIUS*2+2, NODE_RADIUS*2+2);
  return path;
}

void BinaryTreeNode::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *)
{
//  painter->setPen(Qt::NoPen);
//  painter->setBrush(Qt::darkGray);
//  painter->drawEllipse(-3.5, -3.5, 5, 5);

  QRadialGradient gradient(-3, -3, 10);
  gradient.setColorAt(0, m_bHighlighted?QColor(Qt::red).lighter() : Qt::yellow);
  gradient.setColorAt(1, m_bHighlighted?Qt::red : Qt::darkYellow);

  painter->setBrush(gradient);

  painter->setPen(QPen(m_bHighlighted?Qt::darkRed:Qt::black, 0));
  //painter->drawEllipse(-10, -10, 20, 20);
  painter->drawEllipse(-NODE_RADIUS, -NODE_RADIUS, NODE_RADIUS*2, NODE_RADIUS*2);
}

void BinaryTreeNode::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  QGraphicsItem::mousePressEvent(event);
}

void BinaryTreeNode::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
  QGraphicsItem::mouseReleaseEvent(event);
}

void BinaryTreeNode::SetHighlighted(bool bHighlight, bool whole_branch)
{
  m_bHighlighted = bHighlight;
  if (whole_branch)
  {
    foreach (BinaryTreeEdge* edge, edgeList)
    {
      if (edge->destNode() && edge->destNode() != this)
      {
        edge->destNode()->SetHighlighted(bHighlight, whole_branch);
        edge->SetHighlighted(bHighlight);
      }
    }
  }
  update();
}

