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

public slots:
  void OnComboMask(int n);
  void OnLineEditScale(const QString& strg);
  void OnLineEditMaskThresholdChanged();
  void OnLineEditColorThreshold(const QString& strg);
  void OnSliderColorThreshold(int nVal);

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );

private:
  Ui::PanelODF *ui;
  QList<QWidget*> m_colorThresholdWidgets;
};

#endif // PANELODF_H
