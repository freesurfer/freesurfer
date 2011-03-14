#ifndef INFOTREEWIDGET_H
#define INFOTREEWIDGET_H

#include <QTreeWidget>

class QLineEdit;
class QTreeWidgetItem;

class InfoTreeWidget : public QTreeWidget
{
    Q_OBJECT
public:
    InfoTreeWidget(QWidget* parent = 0);

signals:
    void RASChangeTriggered(double x, double y, double z);

protected slots:
    void OnMousePositionChanged();
    void OnCursorPositionChanged();
    void UpdateAll();
    void OnItemClicked(QTreeWidgetItem * item, int column);
    void OnEditFinished();

protected:
    void showEvent(QShowEvent *);
    void keyPressEvent(QKeyEvent *event);
    void mousePressEvent(QMouseEvent *event);

private:
    double m_dRAS[3];
    QLineEdit*  m_editor;
    QTreeWidgetItem* m_itemEdited;
};

#endif // INFOTREEWIDGET_H
