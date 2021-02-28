#ifndef PANELODF_H
#define PANELODF_H

#include "PanelLayer.h"

namespace Ui {
class PanelODF;
}

class PanelODF : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelODF(QWidget *parent = nullptr);
  ~PanelODF();

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );

private:
  Ui::PanelODF *ui;
};

#endif // PANELODF_H
