#ifndef PANELROI_H
#define PANELROI_H

#include "PanelLayer.h"

namespace Ui {
    class PanelROI;
}

class PanelROI : public PanelLayer
{
    Q_OBJECT

public:
    explicit PanelROI(QWidget *parent = 0);
    ~PanelROI();

protected:
    void DoUpdateWidgets();
    void DoIdle();
    virtual void ConnectLayer( Layer* layer );

protected slots:
    void OnSliderOpacity( int val );

private:
    Ui::PanelROI *ui;
};

#endif // PANELROI_H
