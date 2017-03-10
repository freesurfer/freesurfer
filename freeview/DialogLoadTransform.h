#ifndef DIALOGLOADTRANSFORM_H
#define DIALOGLOADTRANSFORM_H

#include <QDialog>

namespace Ui {
class DialogLoadTransform;
}

class DialogLoadTransform : public QDialog
{
  Q_OBJECT

public:
  explicit DialogLoadTransform(QWidget *parent = 0);
  ~DialogLoadTransform();

  QString GetFilename();

  int GetSampleMethod();

public slots:
  void OnOK();
  void OnButtonOpen();

private:
  Ui::DialogLoadTransform *ui;

  QString     m_strLastDir;
};

#endif // DIALOGLOADTRANSFORM_H
