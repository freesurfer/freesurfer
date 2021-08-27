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
#ifndef DIALOGLOADVOLUME_H
#define DIALOGLOADVOLUME_H

#include <QDialog>
#include <QStringList>

namespace Ui
{
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

  bool GetLoadAsVector();

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
