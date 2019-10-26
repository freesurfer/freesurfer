#ifndef WINDOWEDITANNOTATION_H
#define WINDOWEDITANNOTATION_H

#include <QWidget>

#include "colortab.h"

namespace Ui {
class WindowEditAnnotation;
}

class LayerSurface;
class Layer;
class QTreeWidgetItem;

class WindowEditAnnotation : public QWidget
{
  Q_OBJECT

public:
  explicit WindowEditAnnotation(QWidget *parent = nullptr);
  ~WindowEditAnnotation();

  void showEvent(QShowEvent* e);

  int GetCurrentIndex();

signals:
  void LabelClicked(int n);

public slots:
  void OnActiveSurfaceChanged(Layer* layer);
  void UpdateUI();
  void PopulateColorTable(bool bForce = false);
  void OnExistingLabelClicked(QTreeWidgetItem* item);
  void OnExistingLabelItemChanged(QTreeWidgetItem *item);
  void OnCheckBoxShowAllLabels(int);
  void OnSurfaceVertexClicked(LayerSurface* surf);
  void OnEditNameReturnPressed();
  void OnColorChanged(const QColor& color);
  void PopulateAvailableColorTable();
  void OnAvailableLabelClicked(QTreeWidgetItem* item);
  void OnButtonSet();

private:
  void UpdateLabelItem(QTreeWidgetItem* item, int i, const QString& name, const QColor& color);

  Ui::WindowEditAnnotation *ui;

  LayerSurface* m_layerSurface;
  COLOR_TABLE*  m_ctab;
};

#endif // WINDOWEDITANNOTATION_H
