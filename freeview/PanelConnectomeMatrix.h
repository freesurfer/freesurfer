#ifndef PANELCONNECTOMEMATRIX_H
#define PANELCONNECTOMEMATRIX_H

#include "PanelLayer.h"

namespace Ui {
    class PanelConnectomeMatrix;
}

class PanelConnectomeMatrix : public PanelLayer
{
    Q_OBJECT

public:
    explicit PanelConnectomeMatrix(QWidget *parent = 0);
    ~PanelConnectomeMatrix();

protected:
  void DoIdle();
  void DoUpdateWidgets();
  virtual void ConnectLayer( Layer* layer );

private:
    Ui::PanelConnectomeMatrix *ui;
};

#endif // PANELCONNECTOMEMATRIX_H
