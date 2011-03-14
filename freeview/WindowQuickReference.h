#ifndef WINDOWQUICKREFERENCE_H
#define WINDOWQUICKREFERENCE_H

#include <QWidget>

namespace Ui {
    class WindowQuickReference;
}

class WindowQuickReference : public QWidget
{
    Q_OBJECT

public:
    explicit WindowQuickReference(QWidget *parent = 0);
    ~WindowQuickReference();

private:
    Ui::WindowQuickReference *ui;
};

#endif // WINDOWQUICKREFERENCE_H
