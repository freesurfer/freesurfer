#ifndef DIALOGLOADSURFACE_H
#define DIALOGLOADSURFACE_H

#include <QDialog>
#include <QStringList>

namespace Ui {
class DialogLoadSurface;
}

class DialogLoadSurface : public QDialog
{
  Q_OBJECT

public:
  explicit DialogLoadSurface(QWidget *parent = 0);
  ~DialogLoadSurface();

  QString GetFilename();
  QStringList GetSupFiles();

public slots:
  void accept();
  void OnOpen();

private:
  void UpdateStatus();

  Ui::DialogLoadSurface *ui;
};

#endif // DIALOGLOADSURFACE_H
