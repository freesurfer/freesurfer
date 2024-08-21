#ifndef DIALOGLOADODF_H
#define DIALOGLOADODF_H

#include <QDialog>

namespace Ui {
class DialogLoadODF;
}

class DialogLoadODF : public QDialog
{
  Q_OBJECT

public:
  explicit DialogLoadODF(QWidget *parent = nullptr);
  ~DialogLoadODF();

  QString GetOdfFile();
  QString GetVertexFile();
  QString GetMeshFile();

public slots:
  void OnButtonODF();
  void OnButtonVertex();
  void OnButtonMesh();
  void OnCheckboxDTK(bool b);
  void OnOK();

private:
  Ui::DialogLoadODF *ui;
  QString m_strLastDir;
};

#endif // DIALOGLOADODF_H
