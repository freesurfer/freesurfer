#ifndef DIALOGLOADDTI_H
#define DIALOGLOADDTI_H

#include <QDialog>

namespace Ui {
    class DialogLoadDTI;
}

class DialogLoadDTI : public QDialog
{
    Q_OBJECT

public:
    explicit DialogLoadDTI(QWidget *parent = 0);
    ~DialogLoadDTI();

    QString GetVectorFileName();
     QString GetFAFileName();
     QString GetRegFileName();

     bool IsToResample();

     void SetLastDir( const QString& dir )
     {
       m_strLastDir = dir;
     }

     void Initialize( bool bResample, bool bEnableCheckBox );

 protected slots:
     void OnOK();
     void OnButtonVector();
     void OnButtonFA();
     void OnButtonRegistration();

private:
    Ui::DialogLoadDTI *ui;
    QString m_strLastDir;
};

#endif // DIALOGLOADDTI_H
