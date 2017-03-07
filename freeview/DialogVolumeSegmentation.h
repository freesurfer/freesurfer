#ifndef DIALOGVOLUMESEGMENTATION_H
#define DIALOGVOLUMESEGMENTATION_H

#include <QDialog>

namespace Ui {
class DialogVolumeSegmentation;
}

class DialogVolumeSegmentation : public QDialog
{
  Q_OBJECT

public:
  explicit DialogVolumeSegmentation(QWidget *parent = 0);
  ~DialogVolumeSegmentation();


  bool ValidateInput();

public slots:
  void OnButtonRun();
  void OnButtonRestore();
  void UpdateVolumes();

private:
  Ui::DialogVolumeSegmentation *ui;
};

#endif // DIALOGVOLUMESEGMENTATION_H
