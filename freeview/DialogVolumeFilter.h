#ifndef DIALOGVOLUMEFILTER_H
#define DIALOGVOLUMEFILTER_H

#include <QDialog>

namespace Ui {
    class DialogVolumeFilter;
}

class VolumeFilter;

class DialogVolumeFilter : public QDialog
{
    Q_OBJECT

public:
    explicit DialogVolumeFilter(QWidget *parent = 0);
    ~DialogVolumeFilter();

    void SetFilter( VolumeFilter* filter );

    int GetKernelSize();

    double GetSigma();

    void SetSigma( double dvalue );

    void ShowSigma( bool bShow );

protected slots:
    void OnOK();

private:
    Ui::DialogVolumeFilter *ui;
    VolumeFilter* m_filter;
};

#endif // DIALOGVOLUMEFILTER_H
