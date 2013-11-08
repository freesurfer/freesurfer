#ifndef PANELALLLAYERS_H
#define PANELALLLAYERS_H

#include <QScrollArea>
#include <QList>

namespace Ui {
    class PanelAllLayers;
}

class Layer;
class QTreeWidgetItem;
class PanelLayer;

class PanelAllLayers : public QScrollArea
{
  Q_OBJECT

public:
  explicit PanelAllLayers(QWidget *parent = 0);
  ~PanelAllLayers();

  QString GetCurrentLayerType();

  QList<Layer*> GetSelectedLayers(const QString& layerType);

signals:
  void LayerTypeTriggered(const QString& type);

public slots:
  void OnActiveLayerChanged(Layer* curLayer);
  void OnLayerRemoved(Layer* added_layer);
  void OnLayerAdded(Layer* removed_layer);
  void RefreshLayerList(Layer* curLayer = NULL);
  void OnCurrentItemChanged(QTreeWidgetItem* item);
  void OnItemDoubleClicked(QTreeWidgetItem* item);
  void OnItemChanged(QTreeWidgetItem* item);
  void OnItemSelectionChanged();
  void OnLayerChanged();

  void UpdateWidgets();

private:
  void AddLayers(QList<Layer*> layers, const QString& cat_name, Layer* activeLayer, Layer* curLayer = NULL);
  PanelLayer* SetCurrentPanel(const QString& layerType);

  Ui::PanelAllLayers *ui;
};

#endif // PANELALLLAYERS_H
