#ifndef PANELFCD_H
#define PANELFCD_H

#include "PanelLayer.h"

namespace Ui {
    class PanelFCD;
}

class LayerFCD;

class PanelFCD : public PanelLayer
{
    Q_OBJECT

public:
    explicit PanelFCD(QWidget *parent = 0);
    ~PanelFCD();

public slots:
  void OnSliderThresholdReleased();
  void OnSliderSigmaReleased();
  void OnTextThresholdReturned();
  void OnTextSigmaReturned();
  void OnSliderMinAreaReleased();
  void OnTextMinAreaReturned();
  void UpdateLabelList(LayerFCD* layer = NULL);

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );

private:
    Ui::PanelFCD *ui;
};

#endif // PANELFCD_H
