#ifndef WINDOWLAYERINFO_H
#define WINDOWLAYERINFO_H

#include <QWidget>

namespace Ui {
class WindowLayerInfo;
}

class Layer;

class WindowLayerInfo : public QWidget
{
  Q_OBJECT

public:
  explicit WindowLayerInfo(QWidget *parent = nullptr);
  ~WindowLayerInfo();

public slots:
  void UpdateInfo(Layer* layer);

private:
  void Clear();
  void AddLine(const QString& name, const QString& value, bool word_wrap = false);
  void SetCaption(const QString& text);

  Ui::WindowLayerInfo *ui;
};

#endif // WINDOWLAYERINFO_H
