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

  QStringList GetFilenames();
  bool GetLoadAll();

public slots:
  void accept();
  void OnOpen();
  void SetRecentFiles( const QStringList& fns);

private:
  void UpdateStatus();

  Ui::DialogLoadSurface *ui;
};

#endif // DIALOGLOADSURFACE_H
