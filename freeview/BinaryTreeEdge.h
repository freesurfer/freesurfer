#ifndef BINARYTREEEDGE_H
#define BINARYTREEEDGE_H

#include <QGraphicsItem>

class BinaryTreeNode;

class BinaryTreeEdge : public QGraphicsItem
{
public:
    BinaryTreeEdge(BinaryTreeNode *sourceNode, BinaryTreeNode *destNode);

    BinaryTreeNode *sourceNode() const;
    BinaryTreeNode *destNode() const;

    void adjust();

    enum { Type = UserType + 2 };
    int type() const { return Type; }

    void SetHighlighted(bool bHighlight);

protected:
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

private:
    BinaryTreeNode *source, *dest;

    QPointF sourcePoint;
    QPointF destPoint;
    qreal arrowSize;
    bool  m_bHighlighted;
};

#endif // BINARYTREEEDGE_H
