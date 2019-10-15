#ifndef DIALOGNEWANNOTATION_H
#define DIALOGNEWANNOTATION_H

#include <QDialog>

namespace Ui {
class DialogNewAnnotation;
}

class DialogNewAnnotation : public QDialog
{
  Q_OBJECT

public:
  explicit DialogNewAnnotation(QWidget *parent = nullptr, const QString& dir = "");
  ~DialogNewAnnotation();

  QString GetName();
  QString GetColorTableFile();

public slots:
  void OnOK();
  void OnOpen();

private:
  Ui::DialogNewAnnotation *ui;
  QString   m_strDir;
};

#endif // DIALOGNEWANNOTATION_H
