#ifndef PANELTRACK_H
#define PANELTRACK_H

#include "PanelLayer.h"

namespace Ui {
    class PanelTrack;
}

class PanelTrack : public PanelLayer
{
    Q_OBJECT

public:
    explicit PanelTrack(QWidget *parent = 0);
    ~PanelTrack();

protected:
    void DoUpdateWidgets();
    void DoIdle();
    virtual void ConnectLayer( Layer* layer );

private:
    Ui::PanelTrack *ui;
};

#endif // PANELTRACK_H
