#ifndef DIALOGLABELSTATS_H
#define DIALOGLABELSTATS_H

#include <QWidget>

namespace Ui {
    class DialogLabelStats;
}

class DialogLabelStats : public QWidget
{
    Q_OBJECT

public:
    explicit DialogLabelStats(QWidget *parent = 0);
    ~DialogLabelStats();

    void showEvent(QShowEvent *);

  public slots:
    void OnSlicePositionChanged();

private:
    Ui::DialogLabelStats *ui;
};

#endif // DIALOGLABELSTATS_H
