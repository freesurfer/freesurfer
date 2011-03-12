#ifndef DIALOGSAVEPOINTSET_H
#define DIALOGSAVEPOINTSET_H

#include <QDialog>

namespace Ui {
    class DialogSavePointSet;
}

class DialogSavePointSet : public QDialog
{
    Q_OBJECT

public:
    explicit DialogSavePointSet(QWidget *parent = 0);
    ~DialogSavePointSet();

    void SetFileName(const QString& fn);

    void SetLastDir(const QString& dir)
    {
        m_strLastDir = dir;
    }

    void SetType( int nType );

    int GetType();

    QString GetFileName();

protected slots:
    void OnOK();
    void OnOpen();

private:
    Ui::DialogSavePointSet *ui;
    QString m_strLastDir;
};

#endif // DIALOGSAVEPOINTSET_H
