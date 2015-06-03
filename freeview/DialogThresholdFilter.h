#ifndef DIALOGTHRESHOLDFILTER_H
#define DIALOGTHRESHOLDFILTER_H

#include <QDialog>

namespace Ui {
    class DialogThresholdFilter;
}

class DialogThresholdFilter : public QDialog
{
    Q_OBJECT

public:
    explicit DialogThresholdFilter(QWidget *parent = 0);
    ~DialogThresholdFilter();

    void GetThreshold(double* th);

    double GetInValue();

    double GetOutValue();

    bool GetReplaceIn();

    bool GetReplaceOut();

private:
    Ui::DialogThresholdFilter *ui;
};

#endif // DIALOGTHRESHOLDFILTER_H
