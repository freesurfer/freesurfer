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
  friend class PanelSurface;
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
  void UpdateUI(int nIndex = -1);
  void PopulateColorTable(int nIndex = -1);
  void OnExistingLabelClicked(QTreeWidgetItem* item);
  void OnExistingLabelItemChanged(QTreeWidgetItem *item);
  void OnCheckBoxShowAllLabels(int);
  void OnSurfaceVertexClicked(LayerSurface* surf);
  void OnEditColorTextChanged();
  void OnColorChanged(const QColor& color);
  void PopulateAvailableColorTable(bool bForce = false);
  void OnButtonSet();
  void OnButtonFromCTab();
  void OnButtonUndo();
  void OnButtonRedo();
  void OnButtonDelete();
  void OnButtonCleanUp();
  void UpdateActions();
  void OnButtonLoadSegmentation();
  void OnButtonLoadColorTable();

private:
  void UpdateLabelItem(QTreeWidgetItem* item, int i, const QString& name, const QColor& color);
  void UpdateInfoFromItem(QTreeWidgetItem* item);

  Ui::WindowEditAnnotation *ui;

  LayerSurface* m_layerSurface;
  COLOR_TABLE*  m_ctab;
};

#endif // WINDOWEDITANNOTATION_H
