#ifndef DIALOGSAVESCREENSHOT_H
#define DIALOGSAVESCREENSHOT_H

#include <QDialog>
#include "CommonDataStruct.h"

namespace Ui {
    class DialogSaveScreenshot;
}

class DialogSaveScreenshot : public QDialog
{
    Q_OBJECT

public:
    explicit DialogSaveScreenshot(QWidget *parent = 0);
    ~DialogSaveScreenshot();

    QString GetFileName();

    void SetSettings( SettingsScreenshot s );

    SettingsScreenshot GetSettings();

    void SetLastDir( const QString& dir )
    {
      m_strLastDir = dir;
    }

protected slots:
    void OnSave();
    void OnOpen();

private:
    Ui::DialogSaveScreenshot *ui;
    QString m_strLastDir;
};

#endif // DIALOGSAVESCREENSHOT_H
