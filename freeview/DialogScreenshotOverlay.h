#ifndef DIALOGSCREENSHOTOVERLAY_H
#define DIALOGSCREENSHOTOVERLAY_H

#include <QDialog>

namespace Ui {
class DialogScreenshotOverlay;
}

class DialogScreenshotOverlay : public QDialog
{
  Q_OBJECT

public:
  explicit DialogScreenshotOverlay(QWidget *parent = 0);
  ~DialogScreenshotOverlay();

  void showEvent(QShowEvent* e);

public slots:
  void OnButtonOpen();
  void OnButtonSave();
  void OnTimeOut();

private:
  Ui::DialogScreenshotOverlay *ui;

  int   m_nIndex;
};

#endif // DIALOGSCREENSHOTOVERLAY_H
