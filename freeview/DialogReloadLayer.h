#ifndef DIALOGRELOADLAYER_H
#define DIALOGRELOADLAYER_H

#include <QDialog>

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

    int Execute(const QString& layer_name, const QString& layer_type, const QString& fn);
    bool GetCloseLayerFirst();

private:
    Ui::DialogReloadLayer *ui;
};

#endif // DIALOGRELOADLAYER_H
