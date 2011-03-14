#ifndef DIALOGLOADVOLUME_H
#define DIALOGLOADVOLUME_H

#include <QDialog>
#include <QStringList>

namespace Ui {
    class DialogLoadVolume;
}

class DialogLoadVolume : public QDialog
{
    Q_OBJECT

public:
    explicit DialogLoadVolume(QWidget *parent = 0);
    ~DialogLoadVolume();

    QStringList GetVolumeFileNames();

    QString GetRegFileName();

    bool IsToResample();

    int GetSampleMethod();

    QString GetColorMap();

    QString GetLUT();

    void SetLastDir( const QString& dir )
    {
        m_strLastDir = dir;
    }

    void SetRecentFiles( const QStringList& filenames );

protected slots:
    void OnOpen();
    void OnOpenRegistration();
    void OnColorMap( int nSel );
    void OnLUT( int nSel );
    void OnOK();

private:
    void UpdateLUT();

    Ui::DialogLoadVolume *ui;

    QString m_strLastDir;
};

#endif // DIALOGLOADVOLUME_H
