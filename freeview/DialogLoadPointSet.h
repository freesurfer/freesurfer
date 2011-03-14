#ifndef DIALOGLOADPOINTSET_H
#define DIALOGLOADPOINTSET_H

#include <QDialog>

namespace Ui {
    class DialogLoadPointSet;
}

class DialogLoadPointSet : public QDialog
{
    Q_OBJECT
public:
    explicit DialogLoadPointSet(QWidget *parent = 0);
    ~DialogLoadPointSet();

    QStringList GetFileNames();

    int GetPointSetType();

    void SetLastDir( const QString& dir )
    {
      m_strLastDir = dir;
    }

protected slots:
    void OnOK();
    void OnButtonOpen();

private:
    Ui::DialogLoadPointSet *ui;
    QString m_strLastDir;
};

#endif // DIALOGLOADPOINTSET_H
