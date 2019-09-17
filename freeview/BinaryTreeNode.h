#ifndef BINARYTREENODE_H
#define BINARYTREENODE_H


#include <QGraphicsItem>
#include <QList>
#include <QGraphicsView>
#include <QGraphicsPixmapItem>

class BinaryTreeEdge;
class BinaryTreeView;
class QGraphicsSceneMouseEvent;

class BinaryTreeNode : public QGraphicsItem
{
public:
    BinaryTreeNode(BinaryTreeView *graphWidget);

    void addEdge(BinaryTreeEdge *edge);
    QList<BinaryTreeEdge *> edges() const;

    enum { Type = UserType + 1 };
    int type() const { return Type; }

   // void calculateForces();
    bool advance();

    QRectF boundingRect() const;
    QPainterPath shape() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

    void SetHighlighted(bool bHighlight, bool whole_branch = false);

protected:
   // QVariant itemChange(GraphicsItemChange change, const QVariant &value) Q_DECL_OVERRIDE;

    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

private:
    QList<BinaryTreeEdge *> edgeList;
    QPointF newPos;
    bool  m_bHighlighted;
};

#endif // BINARYTREENODE_H
