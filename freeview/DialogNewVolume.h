/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#ifndef DIALOGNEWVOLUME_H
#define DIALOGNEWVOLUME_H

#include <QDialog>

namespace Ui
{
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

  int GetVoxelDataOption();

  int GetDataType();

  LayerMRI* GetTemplate();

protected slots:
  void OnOK();
  void OnToggleCopyVoxelData(bool bCopy);
  void OnToggleVoxelDataOption(bool bChecked);
  void OnToggleMask(bool bMask);

private:
  Ui::DialogNewVolume *ui;
};

#endif // DIALOGNEWVOLUME_H
