#ifndef BUSYINDICATOR_H
#define BUSYINDICATOR_H

#include <QTimer>
#include <QWidget>

class BusyIndicator : public QWidget
{
    Q_OBJECT
public:
    enum IndicatorStyle{StyleRect, StyleEllipse, StyleArc};

    explicit BusyIndicator(QWidget *parent = 0);

    void showEvent(QShowEvent* e);
    void hideEvent(QHideEvent* e);
    void paintEvent(QPaintEvent *);
    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    void setIndicatorStyle(IndicatorStyle);
    void setColor(QColor color);
    IndicatorStyle indicatorStyle() const;

signals:

private slots:
    void rotate();

private:
    QPixmap generatePixmap(int sideLength);
    void drawRectStyle(QPainter *painter);
    void drawEllipseStyle(QPainter *painter);
    void drawArcStyle(QPainter *painter);

    QTimer timer;
    int startAngle;

    IndicatorStyle m_style;

    QColor fillColor;
};

#endif // BUSYINDICATOR_H
