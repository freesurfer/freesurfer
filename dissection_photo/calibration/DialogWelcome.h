#ifndef DIALOGWELCOME_H
#define DIALOGWELCOME_H

#include <QDialog>

namespace Ui {
class DialogWelcome;
}

class DialogWelcome : public QDialog
{
  Q_OBJECT

public:
  explicit DialogWelcome(QWidget *parent = nullptr);
  ~DialogWelcome();

private:
  Ui::DialogWelcome *ui;
};

#endif // DIALOGWELCOME_H
