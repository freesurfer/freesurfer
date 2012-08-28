#ifndef DIALOGLOADSURFACEOVERLAY_H
#define DIALOGLOADSURFACEOVERLAY_H

#include <QDialog>

namespace Ui {
    class DialogLoadSurfaceOverlay;
}

class DialogLoadSurfaceOverlay : public QDialog
{
    Q_OBJECT

public:
    explicit DialogLoadSurfaceOverlay(QWidget *parent = 0);
    ~DialogLoadSurfaceOverlay();

    QString GetFileName();
    QString GetRegistration();

    void SetLastDir(const QString& dir)
    {
      m_strLastDir = dir;
    }

  protected slots:
    void OnOK();
    void OnButtonOpen();
    void OnButtonRegistration();

private:
    Ui::DialogLoadSurfaceOverlay *ui;
    QString m_strLastDir;
};

#endif // DIALOGLOADSURFACEOVERLAY_H
