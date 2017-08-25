#ifndef DIALOGADDPOINTSETSTAT_H
#define DIALOGADDPOINTSETSTAT_H

#include <QDialog>

namespace Ui {
class DialogAddPointSetStat;
}

class DialogAddPointSetStat : public QDialog
{
  Q_OBJECT

public:
  explicit DialogAddPointSetStat(QWidget *parent = 0);
  ~DialogAddPointSetStat();

  QString GetStatName();
  double GetStatValue();

public slots:
  void OnButtonAdd();

private:
  Ui::DialogAddPointSetStat *ui;
};

#endif // DIALOGADDPOINTSETSTAT_H
