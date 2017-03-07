#ifndef DIALOGLOADFCD_H
#define DIALOGLOADFCD_H

#include <QDialog>

namespace Ui {
class DialogLoadFCD;
}

class DialogLoadFCD : public QDialog
{
  Q_OBJECT

public:
  explicit DialogLoadFCD(QWidget *parent = 0);
  ~DialogLoadFCD();

private:
  Ui::DialogLoadFCD *ui;
};

#endif // DIALOGLOADFCD_H
