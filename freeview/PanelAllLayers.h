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
class LayerMRI;

class PanelAllLayers : public QScrollArea
{
  Q_OBJECT

public:
  explicit PanelAllLayers(QWidget *parent = 0);
  ~PanelAllLayers();

  QString GetCurrentLayerType();

  QList<Layer*> GetSelectedLayers(const QString& layerType);

  QList<LayerMRI*> GetLinkedVolumes();

signals:
  void LayerTypeTriggered(const QString& type);
  void ToReorderLayers(const QList<Layer*>& layers);
  void CurrentLayerSelected(Layer* layer);
  void LayerSelectionChanged();

public slots:
  void OnActiveLayerChanged(Layer* curLayer);
  void OnLayerRemoved(Layer* added_layer);
  void OnLayerAdded(Layer* removed_layer);
  void RefreshLayerList(const QList<Layer*>& selectedLayers = QList<Layer*>(), Layer* curLayer = NULL);
  void OnCurrentItemChanged(QTreeWidgetItem* item);
  void OnItemDoubleClicked(QTreeWidgetItem* item);
  void OnItemChanged(QTreeWidgetItem* item);
  void OnItemSelectionChanged();
  void OnLayerChanged();

  void UpdateWidgets();
  void SelectAllLayers();
  void DeselectAllLayers();
  void SetSelectedLayers(const QList<int>& layer_ids);

  PanelLayer* GetPanel(const QString& layer_type);

private:
  void AddLayers(QList<Layer*> layers, const QString& cat_name, Layer* activeLayer,
                 const QList<Layer *>& selectedLayers, Layer* curLayer = NULL);
  PanelLayer* SetCurrentPanel(const QString& layerType);
  void SetItemColor(QTreeWidgetItem* item, double* rgb);

  Ui::PanelAllLayers *ui;
};

#endif // PANELALLLAYERS_H
