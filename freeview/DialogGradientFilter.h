#ifndef DIALOGGRADIENTFILTER_H
#define DIALOGGRADIENTFILTER_H

#include <QDialog>

namespace Ui {
    class DialogGradientFilter;
}

class DialogGradientFilter : public QDialog
{
    Q_OBJECT

public:
    explicit DialogGradientFilter(QWidget *parent = 0);
    ~DialogGradientFilter();

    void SetSmoothing(bool smooth);
    bool GetSmoothing();

    void SetSD(double val);
    double GetSD();

private:
    Ui::DialogGradientFilter *ui;
};

#endif // DIALOGGRADIENTFILTER_H
