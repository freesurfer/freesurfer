#ifndef DIALOGNEWVOLUME_H
#define DIALOGNEWVOLUME_H

#include <QDialog>

namespace Ui {
    class DialogNewVolume;
}

class LayerMRI;

class DialogNewVolume : public QDialog
{
    Q_OBJECT

public:
    explicit DialogNewVolume(QWidget *parent = 0);
    ~DialogNewVolume();

    QString GetVolumeName();
    void SetVolumeName( const QString& name );

    bool GetCopyVoxel();
    void SetCopyVoxel( bool bVoxel );

    int GetDataType();

    LayerMRI* GetTemplate();

protected slots:
    void OnOK();
    void OnToggleCopyVoxelData(bool bCopy);

private:
    Ui::DialogNewVolume *ui;
};

#endif // DIALOGNEWVOLUME_H
