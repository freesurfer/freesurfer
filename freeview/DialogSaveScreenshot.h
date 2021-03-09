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
#ifndef DIALOGSAVESCREENSHOT_H
#define DIALOGSAVESCREENSHOT_H

#include <QDialog>
#include "CommonDataStruct.h"

namespace Ui
{
class DialogSaveScreenshot;
}

class Layer;

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
  void OnSaveLayer();

private:
  Ui::DialogSaveScreenshot *ui;
  QString m_strLastDir;
  QStringList m_listFilenames;
  QList<Layer*> m_listLayers;
  int m_nLayerIndex;
};

#endif // DIALOGSAVESCREENSHOT_H
