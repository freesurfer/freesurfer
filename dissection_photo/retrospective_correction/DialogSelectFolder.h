#ifndef DIALOGSELECTFOLDER_H
#define DIALOGSELECTFOLDER_H

#include <QDialog>

namespace Ui {
class DialogSelectFolder;
}

class DialogSelectFolder : public QDialog
{
  Q_OBJECT

public:
  explicit DialogSelectFolder(QWidget *parent = nullptr);
  ~DialogSelectFolder();

  QString GetInputPath();
  QString GetOutputPath();
  bool IsFourPoint();

public slots:
  void OnButtonInputPath();
  void OnButtonOutputPath();
  void OnButtonRegister();

private:
  Ui::DialogSelectFolder *ui;
};

#endif // DIALOGSELECTFOLDER_H
