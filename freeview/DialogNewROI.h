#ifndef DIALOGNEWROI_H
#define DIALOGNEWROI_H

#include <QDialog>

namespace Ui {
    class DialogNewROI;
}

class LayerMRI;

class DialogNewROI : public QDialog
{
    Q_OBJECT

public:
    explicit DialogNewROI(QWidget *parent = 0);
    ~DialogNewROI();

    QString GetROIName();
    void SetROIName( const QString& name );

    LayerMRI* GetTemplate();

protected slots:
    void OnOK();

private:
    Ui::DialogNewROI *ui;
};

#endif // DIALOGNEWROI_H
