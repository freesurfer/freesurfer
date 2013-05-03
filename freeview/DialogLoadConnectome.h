#ifndef DIALOGLOADCONNECTOME_H
#define DIALOGLOADCONNECTOME_H

#include <QDialog>

namespace Ui {
    class DialogLoadConnectome;
}

class DialogLoadConnectome : public QDialog
{
  Q_OBJECT

public:
  explicit DialogLoadConnectome(QWidget *parent = 0);
  ~DialogLoadConnectome();

  QString GetCMATFilename();
  QString GetParcelFilename();
  QString GetCTABFilename();

public slots:
  void OnButtonOpenCMAT();
  void OnButtonOpenParcel();
  void OnComboBoxColorTable(int n);
  void OnOK();

private:
  void UpdateLUT();

  Ui::DialogLoadConnectome *ui;
  QString m_strLastDir;
};

#endif // DIALOGLOADCONNECTOME_H
