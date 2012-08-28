#ifndef DIALOGSMOOTHSURFACE_H
#define DIALOGSMOOTHSURFACE_H

#include <QDialog>

namespace Ui {
    class DialogSmoothSurface;
}

class DialogSmoothSurface : public QDialog
{
    Q_OBJECT

public:
    explicit DialogSmoothSurface(QWidget *parent = 0);
    ~DialogSmoothSurface();

public slots:
    void OnApply();
    void OnMethod(int nMethod);

private:
    Ui::DialogSmoothSurface *ui;

    bool ValidateAll();
};

#endif // DIALOGSMOOTHSURFACE_H
