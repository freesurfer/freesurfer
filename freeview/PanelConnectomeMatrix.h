#ifndef PANELCONNECTOMEMATRIX_H
#define PANELCONNECTOMEMATRIX_H

#include "PanelLayer.h"

namespace Ui {
    class PanelConnectomeMatrix;
}

class QTreeWidget;

class PanelConnectomeMatrix : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelConnectomeMatrix(QWidget *parent = 0);
  ~PanelConnectomeMatrix();

  void PopulateColorTable();

protected slots:
  void OnCurrentFromChanged();
  void OnCurrentToChanged();

  void OnCheckBoxToAll(bool bChecked);
  void UpdateToLabelVisibility();

  void OnSliderFromOpacity(int val);
  void OnSliderToOpacity(int val);

  void OnLineEditSplineRadius(const QString& strg);

protected:
  void DoIdle();
  void DoUpdateWidgets();
  virtual void ConnectLayer( Layer* layer );

private:
  void AddColorTableItem(int value, const QString& name, const QColor& color, QTreeWidget* treeWidget);

  Ui::PanelConnectomeMatrix *ui;
  bool  m_bColorTableDirty;
};

#endif // PANELCONNECTOMEMATRIX_H
