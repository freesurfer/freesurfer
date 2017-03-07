#ifndef DIALOGSETCAMERA_H
#define DIALOGSETCAMERA_H

#include <QDialog>

namespace Ui {
class DialogSetCamera;
}

class DialogSetCamera : public QDialog
{
  Q_OBJECT

public:
  explicit DialogSetCamera(QWidget *parent = 0);
  ~DialogSetCamera();

public slots:
  void OnReset();
  void OnRefresh();
  void OnSettingChanged();

private:
  Ui::DialogSetCamera *ui;
};

#endif // DIALOGSETCAMERA_H
