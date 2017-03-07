#ifndef DIALOGTHRESHOLDVOLUME_H
#define DIALOGTHRESHOLDVOLUME_H

#include <QDialog>

namespace Ui {
class DialogThresholdVolume;
}

class DialogThresholdVolume : public QDialog
{
  Q_OBJECT

public:
  explicit DialogThresholdVolume(QWidget *parent = 0);
  ~DialogThresholdVolume();

public slots:
  void UpdateVolumes();
  void OnButtonRun();
  void OnButtonReset();

private:
  bool ValidateInputs();

  Ui::DialogThresholdVolume *ui;
};

#endif // DIALOGTHRESHOLDVOLUME_H
