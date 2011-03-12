#ifndef TOOLWINDOWROIEDIT_H
#define TOOLWINDOWROIEDIT_H

#include <QWidget>

namespace Ui {
    class ToolWindowROIEdit;
}

class ToolWindowROIEdit : public QWidget
{
    Q_OBJECT

public:
    explicit ToolWindowROIEdit(QWidget *parent = 0);
    ~ToolWindowROIEdit();

    void UpdateWidgets();

protected:
    virtual void showEvent(QShowEvent *);

protected slots:
    void OnEditMode(QAction *act);

private:
    Ui::ToolWindowROIEdit *ui;
};

#endif // TOOLWINDOWROIEDIT_H
