#ifndef TOOLWINDOWEDIT_H
#define TOOLWINDOWEDIT_H

#include "UIUpdateHelper.h"
#include <QWidget>
#include <QList>

namespace Ui {
    class ToolWindowEdit;
}

class ToolWindowEdit : public QWidget, public UIUpdateHelper
{
    Q_OBJECT

public:
    explicit ToolWindowEdit(QWidget *parent = 0);
    ~ToolWindowEdit();

public slots:
    void UpdateWidgets();

protected slots:
    void OnIdle();
    void OnEditMode( QAction* act );
    void OnComboReference(int sel);
    void OnLineEditContourValue(const QString& strg);
    void OnLineEditSmoothSD(const QString& strg);
    void OnDrawRangeChanged(const QString& strg);
    void OnExcludeRangeChanged(const QString& strg);

protected:
    virtual void showEvent(QShowEvent *);

private:
    Ui::ToolWindowEdit *ui;

    bool m_bToUpdateWidgets;
    QList<QWidget*>  m_widgetsBrushSize;
    QList<QWidget*>  m_widgetsReference;
    QList<QWidget*>  m_widgetsTolerance;
    QList<QWidget*>  m_widgetsConstrain;
    QList<QWidget*>  m_widgetsSmooth;
    QList<QWidget*>  m_widgetsContour;
};

#endif // TOOLWINDOWEDIT_H
