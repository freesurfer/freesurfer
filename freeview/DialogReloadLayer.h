#ifndef DIALOGRELOADLAYER_H
#define DIALOGRELOADLAYER_H

#include <QDialog>
#include <QList>

namespace Ui {
class DialogReloadLayer;
}

class Layer;

class DialogReloadLayer : public QDialog
{
  Q_OBJECT

public:
  explicit DialogReloadLayer(QWidget *parent = 0);
  ~DialogReloadLayer();

  int Execute(const QList<Layer*>& layers);
  bool GetCloseLayerFirst();

private:
  Ui::DialogReloadLayer *ui;
};

#endif // DIALOGRELOADLAYER_H
