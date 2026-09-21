#ifndef DIALOGGIFMAKER_H
#define DIALOGGIFMAKER_H

#include <QDialog>

namespace Ui {
class DialogGifMaker;
}

class GifWriterWrapper;

class DialogGifMaker : public QDialog
{
  Q_OBJECT

public:
  explicit DialogGifMaker(QWidget *parent = nullptr);
  ~DialogGifMaker();

  void hideEvent(QHideEvent* e);

public slots:
  void OnButtonAdd();
  void OnButtonClear()
  {
    Reset();
  }
  void OnButtonSave();
  void OnButtonCopy();
  void OnButtonClose();
  void Reset();

private:
  Ui::DialogGifMaker *ui;
  int m_nNumberOfFrames;

  GifWriterWrapper* m_gif;
  QString  m_strTempFilename;
};

#endif // DIALOGGIFMAKER_H
